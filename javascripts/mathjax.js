window.MathJax = {
    tex: {
        inlineMath: [["\\(", "\\)"]],
        displayMath: [["\\[", "\\]"]],
        processEscapes: true,
        processEnvironments: true,
        packages: {'[+]': ['boldsymbol']}
    },
    options: {
        ignoreHtmlClass: ".*|",
        processHtmlClass: "arithmatex"
    },
    loader: {load: ['[tex]/boldsymbol']}
};

document$.subscribe(() => { 
    MathJax.typesetPromise()
})
