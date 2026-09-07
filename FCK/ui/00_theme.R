# ==============================================================================
# ui/00_theme.R — shared styling
#
# Union of the two source apps' <style> blocks (WaPaa1_3.R lines 63-77 and
# CIRCAREG.R lines 56-61); the class names are referenced from ported tabs of
# both families, so both sets are kept.
# ==============================================================================

ui_theme_css <- tags$head(
  tags$style(HTML("
    .content-wrapper, .right-side {
      background-color: #f4f4f4;
    }
    /* from WaPaa: functional ANOVA / post-hoc tabs */
    .significance-legend {
      background-color: #f0f0f0;
      padding: 10px;
      border-radius: 5px;
      margin: 10px 0;
    }
    .pairwise-summary-box {
      background-color: #e8f4fd;
      padding: 15px;
      border-radius: 8px;
      margin: 10px 0;
      border-left: 4px solid #2196F3;
    }
    /* from CIRCAREG: regression tabs */
    .reg-control-panel {
      background-color: #fff;
      padding: 15px;
      border: 1px solid #ddd;
      border-radius: 5px;
      margin-bottom: 15px;
    }
    .box-header { font-weight: bold; }

    /* ------------------------------------------------------------------ *
     * P15.1: every plot is user-resizable.
     *
     * All 47 plot outputs were emitted with a fixed pixel height, chosen once
     * by whoever wrote that tab. A 300px polar plot or a 24-column heatmap is
     * unreadable at the size that suited the author's screen, and there was no
     * way to change it short of editing the source. The wrapper below is added
     * around each plot at runtime (see the script) and carries the native CSS
     * resize handle, so the size is the reader's choice.
     *
     * max-width keeps a drag from pushing the plot outside its dashboard box;
     * the min-* values stop it being collapsed to nothing by accident.
     * ------------------------------------------------------------------ */
    .fck-resizable {
      resize: both;
      overflow: hidden;
      max-width: 100%;
      min-width: 220px;
      min-height: 140px;
      padding-bottom: 6px;
      position: relative;
    }
    /* the native handle is nearly invisible on a white card, so mark the corner */
    .fck-resizable::after {
      content: '';
      position: absolute;
      right: 2px; bottom: 2px;
      width: 10px; height: 10px;
      border-right: 2px solid #b0bec5;
      border-bottom: 2px solid #b0bec5;
      pointer-events: none;
    }
    .fck-resizable:hover::after { border-color: #607d8b; }
    .fck-resizable > .shiny-plot-output,
    .fck-resizable > .html-widget,
    .fck-resizable > .plotly {
      width: 100% !important;
      height: 100% !important;
    }
  ")),

  tags$script(HTML("
    // P15.1: wrap each plot output in a resizable box and keep the plot filling
    // it. Dependency-free on purpose -- shinyjqui would add a package to an app
    // whose environment is not pinned (see the README on renv).
    (function () {
      function enhance(el) {
        if (!el || el.dataset.fckResizable === '1') return;
        if (el.closest('.fck-resizable')) return;
        var h = el.offsetHeight || parseInt(el.style.height, 10) || 400;
        var box = document.createElement('div');
        box.className = 'fck-resizable';
        box.style.height = h + 'px';
        el.parentNode.insertBefore(box, el);
        box.appendChild(el);
        el.dataset.fckResizable = '1';

        if (typeof ResizeObserver === 'undefined') return;
        var t = null;
        new ResizeObserver(function () {
          // Plotly resizes precisely and immediately.
          if (window.Plotly) {
            var g = el.classList.contains('js-plotly-plot')
                  ? el : el.querySelector('.js-plotly-plot');
            if (g) { try { Plotly.Plots.resize(g); } catch (e) {} }
          }
          // Base-R plotOutput re-renders when Shiny recomputes output sizes,
          // which it does on a window resize. Debounced: the observer fires on
          // every pixel of a drag and a re-render per pixel would be unusable.
          if (t) clearTimeout(t);
          t = setTimeout(function () {
            window.dispatchEvent(new Event('resize'));
          }, 150);
        }).observe(box);
      }

      function scan() {
        document.querySelectorAll(
          '.shiny-plot-output, .html-widget-output, .plotly.html-widget'
        ).forEach(enhance);
      }

      if (document.readyState !== 'loading') scan();
      document.addEventListener('DOMContentLoaded', scan);
      // tabs render lazily, so rescan whenever Shiny paints an output
      if (window.jQuery) {
        jQuery(document).on('shiny:value shiny:visualchange', function () {
          setTimeout(scan, 0);
        });
      }
    })();
  "))
)
