
const FULL_DATA_PATH = window.location.pathname.split("/").slice(0, -1).join("/") + "/"

async function saveCurrentPage() {

  // Clone the current document's HTML
  const clone = document.documentElement.cloneNode(true);

  // Update all input and textarea values
  const inputs = clone.querySelectorAll('input, textarea');
  inputs.forEach(input => {
    if (input.tagName === 'TEXTAREA') {
      input.textContent = input.value;
    } else if (input.type === 'text' || input.type === 'password') {
      input.setAttribute('value', input.value);
    } else if (input.type === 'checkbox' || input.type === 'radio') {
      if (input.checked) {
        input.setAttribute('checked', 'checked');
      } else {
        input.removeAttribute('checked');
      }
    }
  });

  // Create a Blob with the modified HTML
  const doctype = '<!DOCTYPE html>';
  const htmlContent = doctype + '\n' + clone.outerHTML;

  try {
    const handle = await window.showSaveFilePicker({
      suggestedName: 'report.html',
      types: [
        {
          description: 'HTML Files',
          accept: { 'text/html': ['.html'] }
        }
      ]
    });
    const writable = await handle.createWritable();
    await writable.write(htmlContent);
    await writable.close();
    console.log('File saved successfully');
  } catch (err) {
    console.error('Save canceled or failed:', err);
  }
}

function createSaveButton() {
  const saveButton = document.createElement('button');
  saveButton.classList.add('btn', 'btn-primary');
  saveButton.textContent = 'Save report';
  saveButton.style.position = 'fixed';
  saveButton.style.top = '1.5rem';
  saveButton.style.right = '1.5rem';
  saveButton.onclick = saveCurrentPage;
  document.body.appendChild(saveButton);
}

document.addEventListener('DOMContentLoaded', () => {
  createSaveButton();
  document.addEventListener('keydown', (event) => {
    if (event.ctrlKey && event.key === 's') {
      event.preventDefault();
      saveCurrentPage();
    }
  });
});
