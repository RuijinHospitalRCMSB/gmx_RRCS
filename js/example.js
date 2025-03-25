// Generate dynamic sidebar navigation
function generateSidebar() {
    const headings = document.querySelectorAll('h1, h2, h3');
    const sidebar = document.getElementById('sidebarNav');
    
    // 清空现有内容
    sidebar.innerHTML = '';
    
    // 获取当前页面路径
    const currentPath = window.location.pathname;
    const isIndexPage = currentPath.includes('index.html') || currentPath.endsWith('/');
    
    // 如果是索引页面，使用静态导航
    if (isIndexPage) {
        const navItems = [
            { id: 'intro', text: 'Introduction' },
            { id: 'beginner', text: 'Beginner Tutorials' },
            { id: 'intermediate', text: 'Intermediate Tutorials' },
            { id: 'advanced', text: 'Advanced Tutorials' }
        ];
        
        navItems.forEach(item => {
            const li = document.createElement('li');
            const a = document.createElement('a');
            a.href = `#${item.id}`;
            a.textContent = item.text;
            li.appendChild(a);
            sidebar.appendChild(li);
        });
    } else {
        // 对于教程页面，动态生成导航
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

// Initialize
window.addEventListener('DOMContentLoaded', () => {
    generateSidebar();
}); 