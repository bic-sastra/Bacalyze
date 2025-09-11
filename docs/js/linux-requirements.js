const searchField = document.querySelector('input');
const tableCells = document.querySelectorAll('table tbody tr td');

searchField.addEventListener('input', () => {
    const searchTerm = searchField.value.toLowerCase();
    for (const cell of tableCells) {
        const row = cell.parentElement;
        const match = cell.textContent.toLowerCase().includes(searchTerm);
        row.style.display = match ? '' : 'none';
    }
});
