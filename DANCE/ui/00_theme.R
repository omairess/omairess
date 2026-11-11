# ==============================================================================
# ui/00_theme.R — shared styling
#
# Union of the two source apps' <style> blocks (WaPaa1_3.R lines 63-77 and
# CIRCAREG.R lines 56-61); the class names are referenced from ported tabs of
# both families, so both sets are kept.
# ==============================================================================

# the explanation helpers (dance_help, dance_status, and helpText routed to
# dance_help) are defined once, in server/00a_help.R, and loaded here too so
# the UI files use the same ones
local({
  cand <- c(if (exists("dance_dir", mode = "function")) file.path(dance_dir("server"), "00a_help.R"),
            "server/00a_help.R", "DANCE/server/00a_help.R", "../server/00a_help.R")
  f <- cand[file.exists(cand)][1]
  if (!is.na(f)) eval(parse(f, encoding = "UTF-8"), envir = globalenv())
})

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
    //
    // AUDIT (P16). The first version of this was measured, after the user
    // reported the whole GUI had got slower, and it was the cause. Three
    // separate mistakes, all mine:
    //
    //   scan() did a full document.querySelectorAll and ran on EVERY
    //   shiny:value -- an event that fires once per output, for all 159 of
    //   them, on every render. Measured on a page with 40 plots and 60 text
    //   outputs: 301 full scans to load the page, and 100 more per refresh.
    //
    //   The ResizeObserver dispatched a GLOBAL window resize, which makes Shiny
    //   recompute sizes for every bound output and re-render every base-R plot.
    //   So dragging one plot re-rendered all of them.
    //
    //   ResizeObserver fires once when observation starts, so simply creating
    //   the wrappers fired 40 of those global resizes during page load.
    //
    // Now: shiny:value enhances only the element that just rendered; a full
    // scan happens on load and on tab changes, coalesced into one animation
    // frame; a plotly widget is resized directly and never through a global
    // event; and only a box actually containing a base-R plot falls back to the
    // window event, since that is the one output kind Shiny re-measures that
    // way. The first observer callback per box is ignored.
    (function () {
      var SEL = '.shiny-plot-output, .html-widget-output, .plotly.html-widget';

      function enhance(el) {
        if (!el || !el.matches || !el.matches(SEL)) return;
        if (el.dataset.fckResizable === '1') return;
        if (el.closest('.fck-resizable')) return;

        var h = el.offsetHeight || parseInt(el.style.height, 10) || 400;
        var box = document.createElement('div');
        box.className = 'fck-resizable';
        box.style.height = h + 'px';
        el.parentNode.insertBefore(box, el);
        box.appendChild(el);
        el.dataset.fckResizable = '1';

        if (typeof ResizeObserver === 'undefined') return;

        // Does this box hold a base-R plot? Only those need Shiny to re-measure.
        var isRPlot = el.classList.contains('shiny-plot-output');
        var first = true, t = null;

        new ResizeObserver(function () {
          if (first) { first = false; return; }   // the initial observe() call

          // THE WIDGET HAS TO BE TOLD, not just the box around it.
          // Plotly.Plots.resize() measures the graph's OWN element, and every
          // plotlyOutput in this app is emitted with a fixed inline
          // height:NNNpx. So dragging the wrapper grew the box, plotly
          // re-measured the unchanged element, and the drawing stayed its
          // original size inside a larger frame -- which reads as the resize
          // covering part of the graph rather than doing nothing. (No quote
          // marks in here: this JS lives inside an R string, and a stray one
          // ends the string mid-comment.) Copy the box size onto the element
          // first, then resize.
          if (window.Plotly) {
            var g = el.classList.contains('js-plotly-plot')
                  ? el : el.querySelector('.js-plotly-plot');
            if (g) {
              el.style.height = box.clientHeight + 'px';
              el.style.width  = '100%';
              if (g !== el) { g.style.height = box.clientHeight + 'px'; g.style.width = '100%'; }
              try { Plotly.Plots.resize(g); } catch (e) {}
              return;
            }
          }
          if (!isRPlot) return;

          // Debounced: the observer fires on every pixel of a drag, and a
          // re-render per pixel is unusable.
          if (t) clearTimeout(t);
          t = setTimeout(function () {
            window.dispatchEvent(new Event('resize'));
          }, 200);
        }).observe(box);
      }

      var pending = false;
      function scanSoon() {
        if (pending) return;
        pending = true;
        requestAnimationFrame(function () {
          pending = false;
          document.querySelectorAll(SEL).forEach(enhance);
        });
      }

      if (document.readyState !== 'loading') scanSoon();
      document.addEventListener('DOMContentLoaded', scanSoon);
      if (window.jQuery) {
        // one element, not one scan, for the common case
        jQuery(document).on('shiny:value', function (e) {
          if (e.target && e.target.nodeType === 1) enhance(e.target);
        });
        // tabs render lazily; a whole panel appears at once, so scan then
        jQuery(document).on('shiny:visualchange shown.bs.tab', scanSoon);
      }
    })();
  ")),
  # ---- EXPLANATIONS BEHIND AN (i) (server/00a_help.R: dance_help) ------------
  # The (i) moves next to the label of the control just above it when both sit
  # in the same container (so a note inside a conditional panel still shows
  # and hides with it); elsewhere it stays where it is, with 'About this'.
  # Clicking opens the text in a pop-up titled with that label. (No double
  # quotes in the JS: it lives inside an R string.)
  tags$style(HTML("
    .dance-help { display: inline-block; margin: 2px 0 6px 0; font-size: 12px; }
    .dance-help-icon { color: #1f5f7a; text-decoration: none; cursor: pointer; }
    .dance-help-icon:hover { color: #0b3c50; text-decoration: none; }
    .dance-help-more { margin-left: 4px; color: #6b7b85; }
    .dance-help-inline { margin: 0 0 0 6px; font-size: 13px; }
    .dance-help-inline .dance-help-more { display: none; }
    #dance-help-modal .modal-body { font-size: 14px; line-height: 1.5; }
  ")),
  tags$script(HTML("
    (function () {
      function modal() {
        var m = document.getElementById('dance-help-modal');
        if (m) return m;
        m = document.createElement('div');
        m.id = 'dance-help-modal'; m.className = 'modal fade'; m.tabIndex = -1;
        m.innerHTML = '<div class=\\'modal-dialog\\'><div class=\\'modal-content\\'>' +
          '<div class=\\'modal-header\\'><button type=\\'button\\' class=\\'close\\' data-dismiss=\\'modal\\'>&times;</button>' +
          '<h4 class=\\'modal-title\\'></h4></div><div class=\\'modal-body\\'></div></div></div>';
        document.body.appendChild(m);
        return m;
      }
      window.danceHelp = function (a, ev) {
        if (ev) { ev.preventDefault(); ev.stopPropagation(); }
        var w = a.closest('.dance-help'); if (!w) return false;
        var body = w.querySelector('.dance-help-body');
        var title = w.getAttribute('data-title');
        if (!title) {
          var lab = w.parentElement;
          if (lab && w.classList.contains('dance-help-inline')) {
            var c = lab.cloneNode(true);
            c.querySelectorAll('.dance-help').forEach(function (x) { x.remove(); });
            title = c.textContent;
          } else title = 'Explanation';
        }
        var m = modal();
        m.querySelector('.modal-title').textContent = (title || 'Explanation').trim().replace(/[:\\s]+$/, '');
        m.querySelector('.modal-body').innerHTML = body ? body.innerHTML : '';
        if (window.jQuery) jQuery(m).modal('show');
        return false;
      };
      function place(el) {
        if (el.getAttribute('data-placed')) return;
        el.setAttribute('data-placed', '1');
        var prev = el.previousElementSibling, lab = null;
        if (prev) {
          if (prev.matches('label, .control-label, h3, h4, h5')) lab = prev;
          else lab = prev.querySelector('label.control-label, .control-label, .checkbox label, label');
        }
        if (lab && !lab.querySelector('.dance-help')) {
          el.classList.add('dance-help-inline');
          lab.appendChild(el);
        }
      }
      function scan(root) {
        (root || document).querySelectorAll('.dance-help:not([data-placed])').forEach(place);
      }
      function start() {
        scan(document);
        new MutationObserver(function (ms) {
          ms.forEach(function (mu) {
            mu.addedNodes.forEach(function (n) {
              if (n.nodeType !== 1) return;
              if (n.classList && n.classList.contains('dance-help')) place(n);
              if (n.querySelectorAll) scan(n);
            });
          });
        }).observe(document.body, { childList: true, subtree: true });
      }
      if (document.readyState !== 'loading') start();
      else document.addEventListener('DOMContentLoaded', start);
    })();
  "))
)
