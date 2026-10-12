// ============================================================================
// rd-logo.js — tint the ".jl" extension of the sidebar wordmark teal.
//
// Documenter renders the sitename as a bare text node inside the brand link:
//   <div class="docs-package-name"><span class="docs-autofit">
//     <a href="…">ReactiveDynamics.jl</a></span></div>
// There is no per-glyph markup to hook a colour onto, so we split the trailing
// ".jl" into its own <span class="rd-jl"> at load time; rd-theme.css paints that
// span with the house teal accent (echoing the identity board's `.jl` treatment).
//
// Guardrails: only the ".jl" SUFFIX is wrapped, we bail if the structure isn't
// the single-text-node shape we expect, and the wrap is idempotent (skipped once
// an .rd-jl span is present), so a re-run — e.g. Documenter's live reload — is a
// no-op rather than a double-wrap.
// ============================================================================
(function () {
  "use strict";

  function tintExtension() {
    var link = document.querySelector(
      "#documenter .docs-sidebar .docs-package-name a"
    );
    if (!link) return;
    if (link.querySelector(".rd-jl")) return; // already wrapped — idempotent

    // Only touch the simple case: a single text node ending in ".jl". Anything
    // richer (already-marked-up name) is left exactly as Documenter emitted it.
    if (link.childNodes.length !== 1) return;
    var node = link.firstChild;
    if (node.nodeType !== Node.TEXT_NODE) return;

    var text = node.nodeValue;
    var ext = ".jl";
    if (!text.endsWith(ext)) return;

    var stem = document.createTextNode(text.slice(0, -ext.length));
    var jl = document.createElement("span");
    jl.className = "rd-jl";
    jl.textContent = ext;

    link.replaceChild(stem, node);
    link.appendChild(jl);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", tintExtension);
  } else {
    tintExtension();
  }
})();
