## AT Web site

The AT web site is handled by [GitHub Pages](https://pages.github.com).
It is generated from the contents of the `docs` directory by a Web site generator
called [jekyll](https://jekyllrb.com). The site is generated at each commit on
the master branch of the AT repository.

As a consequence it is not possible to view the results of a modification of the
Web documentation before merging it into the master branch. The only way of previewing 
is to install [jekyll](https://jekyllrb.com) locally on your machine and run it on your
working copy.

### Generation of the Web site
[jekyll](https://jekyllrb.com) will ignore the subdirectories of the `docs` directory
with a name starting with `.`, `_`, `#` or `~`, move the other ones to the web site
and translate their contents to html.

The layout of the site is based on predefined template. Some are proposed by GitHub Pages,
other ones are available on-line.

[jekyll](https://jekyllrb.com) accepts either html or markdown as input. The simplest 
way of writing is therefore to use the [jekyll markdown flavour](https://daringfireball.net/projects/markdown/syntax).

### Structure of the AT `docs` directory
The docs directory consists in:

* the `_config.yml` file: it defines the template to be used and some default
parameters of the site,
* `_*` directories: they set customizations of the template specific for the AT site,
* 4 header markdown files for the top-level of the site: `index.md`, `about.md`,
  `matlab.md` and `python.md`,
* the `p` directory: contains python-specific pages,
* the `m` directory: contains Matlab-specific pages,
* the `matlab` directory: it contains the html files for the Matlab help system.
These are not part of the Web site, however they are located here to allow some Web page
to include them.

### The AT template
The AT Web site is based on the [Documentation Theme for jekyll](https://idratherbewriting.com/documentation-theme-jekyll/)
by Tom Johnson, with some modifications. It was chosen for its top navigation bar and its possibility of
optional side bars. There is a side bar for python pages, defined by the YAML file
`docs/_data/sidebars/python_sidebar.yml` and another one for Matlab pages defined in
`docs/_data/sidebars/matlab_sidebar.yml`.

### Writing AT documentation
The easiest way of contributing to the AT documentation is to write a markdown file.
Nice formatting can be obtained by making use of the many specific features of
the [jekyll markdown](https://daringfireball.net/projects/markdown/syntax) syntax.
Formulæ and figures can be easily included. Then put this file into the `docs/p`
directory and add a link to it in the sidebar configuration file.

> Due to some restrictions of the theme, a jekyll-specific section delimited by
> "`---`" lines must be present at the top of the file. It does not appear in the output.
> The `permalink` directive defines the location of the html file in the site, which
> must be at the top-level. To avoid name conflicts, please follow the convention
> of prepending a "p" (or "m") to the markdown file name (see on existing files).

For creating markdown files you can also look at [pandoc](https://pandoc.org) which
can convert from a lot of formats (LaTeX, .docx, .odt…) to markdown.