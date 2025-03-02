function hideBrokenLinks() {
  // Find all elements with class 'hide-broken' and hide any child elements with
  // an href that returns an error code
  const elements = document.getElementsByClassName('hide-broken');
  Array.from(elements).forEach(element => {
    // Get all child elements with an href attribute
    const children = element.querySelectorAll('[href]');
    const hrefElements = [element].concat(Array.from(children));
    hrefElements.forEach(el => {
      // Check if the href returns an error code
      if (!el.href) {
        return;
      }
      hideIfLinkBroken(element, el.href);
    });
  });
}

function hideIfLinkBroken(element, url) {
  console.log("Testing element:" + element);
  console.log("Testing href:" + url);
  let testImg = new Image();
  testImg.src = url;
  testImg.onload = () => {
      // Do nothing, link is valid
  };
  testImg.onerror = () => {
      element.style.display = "none";
      console.log("Hiding element with broken link: " + url);
  };
}

document.addEventListener('DOMContentLoaded', hideBrokenLinks);
