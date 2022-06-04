---
title: Sample page
summary: This demonstrates some features of the AT theme
pict2: assets/images/output_69_0.png
hide_sidebar: true
---
If present, the **title** and **summary** defined in the frontmatter are automatically
displayed.\
The table of contents is automatically generated unless "toc: false" is specified.\
A side-bar is displayed unless "hide_sidebar: true" is specified (as here).

## Markdown features
{: #idfn .red}

**In-line formula:** $\eta_c = 1/\gamma^2 - \alpha_c$

**Display formula:**

$$\gamma=\frac{1+\alpha^2}{\beta}$$

## Links
### Using Liquid tags
Use links relative to the root directory and refer to the html file. You can use in-line or referenced links
The referenced link may be anywhere (here at the bottom of the file).

in-line: [Python installation]({{ "p/Installation.html" | relative_url }})\
referenced: [Matlab installation]

**link tag**: {%link p/Installation.md %}\
The link tag cannot be used on GitHub pages because it still does not prepend
site.baseurl.

### Using markdown syntax
Use relative links (absolute links only work from the root directory), refer to the markdown file.
You can use in-line or referenced links.

in-line: [Matlab installation](m/Installation.md)\
referenced: [Python installation]

## Images
Images should be put in the `/assets/images directory`.
### Using Liquid tags
Use links relative to the root directory. You can use in-line or referenced links.

**Image with in-line link:**

![Figure 1]({{ "assets/images/output_33_1.png" | relative_url }})

**Image with referenced link:**

![Figure 2]

### Using markdown syntax

**Image with in-line link:**

![Figure 3](assets/images/output_33_1.png)

**Image with referenced link:**\
The referenced link may be anywhere (here at the bottom of the file)

![Figure 4]

## Theme features

**Image with format and caption** (specific for the AT theme):
{%include image.html src="assets/images/output_67_1.png"
alt="Figure 3" caption="Figure 3: test caption" %}

**Alerts:**
{% include note.html content="Alerts are displayed in colour boxes" %}

**Comments:** (do not appear)
{% comment %} The following lines are the referenced links
and do not appear in the output{% endcomment %}

## Using variables {#location}

{% assign pict1 = "assets/images/output_33_1.png" %}

**site.title**: {{ site.title }}\
**site.description**: {{ site.description }}\
**site.url**: {{ site.url }}\
**site.baseurl**: {{ site.baseurl }}\
**site.github.owner_name**: {{ site.github.owner_name }}\
**site.github.owner_url**: {{ site.github.owner_url }}\
**site.github.owner_gravatar_url**: {{ site.github.owner_gravatar_url }}\
**site.github.repository_url**: {{ site.github.repository_url }}\
**pict1**: {{ pict1 }}\
**page.pict2**: {{ page.pict2 }}

{% comment %} Link references (do not appear in the output) {% endcomment %}
[Figure 2]: {{ page.pict2 | relative_url }}
[Figure 4]: assets/images/output_69_0.png
[Matlab installation]: {{ "m/Installation.html" | relative_url }}
[Python installation]: p/Installation.md
