'use strict';

const { app, BrowserWindow, ipcMain, screen, dialog, powerSaveBlocker } = require('electron');
const path = require('path');
const fs = require('fs');
const os = require('os');

/* ------------------------------------------------------------------
 * Timing-related switches.
 *
 * We deliberately do NOT disable vsync. Frame-locked presentation is the
 * basis of the whole timing model: the task counts real vsync intervals,
 * so tearing the renderer loose from the display would make onset times
 * less knowable, not more.
 *
 * What we do disable is every mechanism Chromium uses to throttle or
 * deprioritise a window it thinks is idle or hidden.
 * ------------------------------------------------------------------ */
app.commandLine.appendSwitch('disable-background-timer-throttling');
app.commandLine.appendSwitch('disable-renderer-backgrounding');
app.commandLine.appendSwitch('disable-backgrounding-occluded-windows');
app.commandLine.appendSwitch('disable-features', 'CalculateNativeWinOcclusion');
if (process.platform === 'darwin') {
  app.commandLine.appendSwitch('force_high_performance_gpu');
}

let win = null;
let blockerId = null;

function createWindow() {
  win = new BrowserWindow({
    width: 1100,
    height: 820,
    backgroundColor: '#000000',
    show: false,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      nodeIntegration: false,
      // Critical: without this, Chromium throttles rAF and timers whenever
      // the window is occluded or backgrounded, which would silently corrupt
      // a 40-minute trial.
      backgroundThrottling: false
    }
  });

  win.loadFile(path.join(__dirname, 'renderer', 'index.html'));
  win.once('ready-to-show', () => win.show());
  win.on('closed', () => { win = null; });
}

app.whenReady().then(() => {
  createWindow();
  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) createWindow();
  });
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') app.quit();
});

app.on('will-quit', () => {
  if (blockerId !== null && powerSaveBlocker.isStarted(blockerId)) {
    powerSaveBlocker.stop(blockerId);
  }
});

/* ---------------- IPC ---------------- */

ipcMain.handle('display-info', () => {
  const primary = screen.getPrimaryDisplay();
  const current = win
    ? screen.getDisplayMatching(win.getBounds())
    : primary;

  const describe = (d) => ({
    id: d.id,
    bounds: d.bounds,
    scaleFactor: d.scaleFactor,
    // Electron reports the OS-declared refresh rate. It is the *nominal*
    // rate, and can disagree with what the compositor actually achieves,
    // which is exactly why the renderer measures it independently.
    displayFrequency: typeof d.displayFrequency === 'number' ? d.displayFrequency : null,
    colorDepth: d.colorDepth,
    internal: d.internal,
    label: d.label || null
  });

  return {
    current: describe(current),
    isPrimary: current.id === primary.id,
    displayCount: screen.getAllDisplays().length,
    platform: process.platform,
    arch: process.arch,
    osRelease: os.release(),
    cpuModel: (os.cpus()[0] || {}).model || null,
    electron: process.versions.electron,
    chrome: process.versions.chrome,
    appVersion: app.getVersion()
  };
});

ipcMain.handle('set-fullscreen', (_e, on) => {
  if (!win) return false;
  win.setFullScreen(!!on);
  return win.isFullScreen();
});

ipcMain.handle('power-block', (_e, on) => {
  if (on) {
    if (blockerId === null || !powerSaveBlocker.isStarted(blockerId)) {
      blockerId = powerSaveBlocker.start('prevent-display-sleep');
    }
    return true;
  }
  if (blockerId !== null && powerSaveBlocker.isStarted(blockerId)) {
    powerSaveBlocker.stop(blockerId);
  }
  blockerId = null;
  return false;
});

ipcMain.handle('save-csv', async (_e, suggestedName, contents) => {
  if (!win) return { ok: false, reason: 'no window' };
  const res = await dialog.showSaveDialog(win, {
    title: 'Save data',
    defaultPath: path.join(app.getPath('documents'), suggestedName),
    filters: [{ name: 'CSV', extensions: ['csv'] }]
  });
  if (res.canceled || !res.filePath) return { ok: false, reason: 'canceled' };
  try {
    fs.writeFileSync(res.filePath, contents, 'utf8');
    return { ok: true, path: res.filePath };
  } catch (err) {
    return { ok: false, reason: err.message };
  }
});
