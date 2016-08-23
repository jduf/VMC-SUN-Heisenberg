# rst2html-style

This is a css stylesheet for making better looking html output from [reStructuredText](http://docutils.sourceforge.net/rst.html) and [rst2html](http://docutils.sourceforge.net/docs/user/tools.html#rst2html-py).

Take a look at [this sample](http://cskeeters.bitbucket.org/rst2html-style/doc.html).

# Pygments Hack

In order for code coloring to work, the rst2html (or rst2html.py) executable must be modified to insert pygements classes into the html output.  This code is documented [here](http://stackoverflow.com/a/14880062/319894).

# Generating html

  rst2html --stylesheet=style.css doc.txt > doc.html

