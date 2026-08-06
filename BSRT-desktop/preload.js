'use strict';

const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('bsrt', {
  getDisplayInfo: () => ipcRenderer.invoke('display-info'),
  setFullscreen: (on) => ipcRenderer.invoke('set-fullscreen', on),
  preventDisplaySleep: (on) => ipcRenderer.invoke('power-block', on),
  saveCsv: (name, contents) => ipcRenderer.invoke('save-csv', name, contents)
});
