// Generate dynamic sidebar navigation
function generateSidebar() {
    const headings = document.querySelectorAll('h1, h2, h3');
    const sidebar = document.getElementById('sidebarNav');
    
    // 清空现有内容
    sidebar.innerHTML = '';
    
    headings.forEach(heading => {
        if (!heading.id) return;
        
        const listItem = document.createElement('li');
        const link = document.createElement('a');
        link.href = `#${heading.id}`;
        link.textContent = heading.textContent;
        
        // 添加基于标题级别的缩进
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
        if (target) {
            target.scrollIntoView({ 
                behavior: 'smooth',
                block: 'start'
            });
        }
    });
});

// 代码高亮函数
function initializeCodeHighlighting() {
    if (typeof hljs !== 'undefined') {
        document.querySelectorAll('pre code').forEach(block => {
            // 检测代码内容确定语言
            let language = 'bash';
            const content = block.textContent || '';
            
            if (content.includes('import') || 
                content.includes('params =') || 
                content.includes('def ')) {
                language = 'python';
            }
            
            // 设置语言类
            block.classList.add(`language-${language}`);
            block.parentElement.setAttribute('data-language', language);
            
            // 应用高亮
            hljs.highlightElement(block);
        });
        console.log('Code highlighting applied');
    } else {
        console.warn('highlight.js not loaded');
    }
}

// Initialize
window.addEventListener('DOMContentLoaded', () => {
    generateSidebar();
    
    // 等待一小段时间确保 highlight.js 已加载
    setTimeout(initializeCodeHighlighting, 100);
});

// Load the navbar
document.addEventListener('DOMContentLoaded', function() {
    fetch('../Tutorials/navbar.html')
        .then(response => response.text())
        .then(data => {
            document.querySelector('body').insertAdjacentHTML('afterbegin', data);
        })
        .catch(error => console.error('Error loading navbar:', error));
});