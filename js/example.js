// Generate dynamic sidebar navigation
function generateSidebar() {
    const headings = document.querySelectorAll('h1, h2, h3');
    const sidebar = document.getElementById('sidebarNav');
    
    headings.forEach(heading => {
        if (!heading.id) return;
        
        const listItem = document.createElement('li');
        const link = document.createElement('a');
        link.href = `#${heading.id}`;
        link.textContent = heading.textContent;
        
        // Add indentation based on heading level
        switch(heading.tagName) {
            case 'H2':
                listItem.style.paddingLeft = '20px';
                break;
            case 'H3':
                listItem.style.paddingLeft = '40px';
                break;
        }
        
        listItem.appendChild(link);
        sidebar.appendChild(listItem);
    });
}

// Smooth scrolling
document.querySelectorAll('a[href^="#"]').forEach(anchor => {
    anchor.addEventListener('click', function(e) {
        e.preventDefault();
        const target = document.querySelector(this.hash);
        target.scrollIntoView({ 
            behavior: 'smooth',
            block: 'start'
        });
    });
});

// Initialize
window.addEventListener('DOMContentLoaded', () => {
    generateSidebar();
}); 