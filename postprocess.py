import re

b = open('surfacewaves.html')
bs = b.read()
b.close()
bs = re.sub('<pre>', '<pre data-executable="true" data-language="python">', bs)

s = """
<script type="text/x-thebe-config">
  {
    requestKernel: true,
    bootstrap: false,
    binderOptions: {
      repo: "mikaem/surface-waves-on-liquid-films",
      ref: "master",
      binderUrl: "https://mybinder.org",
      // select repository source (optional). Supports Github(default), Gitlab, and Git
      repoProvider: "github",
    },
  }
</script>

<script src="https://unpkg.com/thebelab@0.4.0/lib/index.js"></script>

<script type="text/x-mathjax-config">
    MathJax.Hub.Config({
        tex2jax: {
            inlineMath: [ ['$','$'], ["\\(","\\)"] ],
            displayMath: [ ['$$','$$'], ["\\[","\\]"] ],
            processEscapes: true,
            processEnvironments: true
        },
        // Center justify equations in code and markdown cells. Elsewhere
        // we use CSS to left justify single line equations in code cells.
        displayAlign: 'center',
        "HTML-CSS": {
            styles: {'.MathJax_Display': {"margin": 0}},
            linebreaks: { automatic: true }
        },
        TeX: {
          equationNumbers: { autoNumber: "AMS" , useLabelIds: true },
          extensions: ["cancel.js", "AMSmath.js"]
        }
    });
</script>

<button id="activateButton"  style="width: 150px; height: 75px; font-size: 1.5em;">Activate</button>
<script>
var bootstrapThebe = function() {
    thebelab.bootstrap();
}

document.querySelector("#activateButton").addEventListener('click', bootstrapThebe)
</script>

"""
b = open('surfacewaves.html', 'w')
b.write(bs[:270]+s+bs[270:])
b.close()
