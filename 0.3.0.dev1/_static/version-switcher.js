(function () {
  "use strict";

  function metaContent(name, fallback) {
    const element = document.querySelector(`meta[name="${name}"]`);
    return element ? element.content : fallback;
  }

  function addSwitcher() {
    const container = document.querySelector(".wy-side-nav-search");
    if (!container || container.querySelector(".edrixs-version-switcher")) {
      return;
    }

    const builtVersion = metaContent("edrixs-doc-version", "unknown");
    const siteRoot = metaContent("edrixs-doc-site-root", "/edrixs/");
    const rootUrl = new URL(siteRoot, window.location.origin);
    const wrapper = document.createElement("div");
    wrapper.className = "edrixs-version-switcher";

    const label = document.createElement("label");
    label.htmlFor = "edrixs-doc-version-select";
    label.textContent = `Version: ${builtVersion}`;
    wrapper.appendChild(label);
    container.appendChild(wrapper);

    fetch(new URL("versions.json", rootUrl))
      .then(function (response) {
        if (!response.ok) {
          throw new Error(`Unable to load versions: ${response.status}`);
        }
        return response.json();
      })
      .then(function (manifest) {
        if (!Array.isArray(manifest.versions) || manifest.versions.length === 0) {
          return;
        }

        const relativePath = window.location.pathname.startsWith(rootUrl.pathname)
          ? window.location.pathname.slice(rootUrl.pathname.length)
          : "";
        const pathParts = relativePath.split("/").filter(Boolean);
        const publishedNames = manifest.versions.map(function (item) {
          return item.name;
        });
        const pathVersion = publishedNames.includes(pathParts[0])
          ? pathParts.shift()
          : builtVersion;
        const pagePath = pathParts.join("/");

        const select = document.createElement("select");
        select.id = "edrixs-doc-version-select";
        select.setAttribute("aria-label", "Documentation version");

        manifest.versions.forEach(function (item) {
          const option = document.createElement("option");
          option.textContent = item.title;
          option.value = new URL(`${item.url}${pagePath}`, rootUrl).href;
          option.selected = item.name === pathVersion;
          select.appendChild(option);
        });

        select.addEventListener("change", function () {
          window.location.assign(select.value);
        });
        label.textContent = "Documentation version";
        wrapper.appendChild(select);
      })
      .catch(function () {
        // A standalone local build has no site-wide versions manifest. The
        // package version label above remains useful in that case.
      });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", addSwitcher);
  } else {
    addSwitcher();
  }
})();
