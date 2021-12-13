---
title: Sample page
summary: This demonstrates some features of the AT theme
pict2: assets/images/output_69_0.png
---
{% assign pict1 = "assets/images/output_69_0.png" %}

If present, the **title** and **summary** defined in the frontmatter are automatically
displayed.\
The table of contents is automatically generated unless "toc: false" is specified.

## Markdown features
{: #idfn .red}

**In-line formula:** $$\eta_c = 1/\gamma^2 - \alpha_c$$

**Display formula:**

$$\gamma=\frac{1+\alpha^2}{\beta}$$

**link:**
[Python installation]({{ "p/Installation.html" | relative_url }})

**Image with in-line link:**

![Figure 1]({{ "assets/images/output_33_1.png" | relative_url }})

**Image with referenced link:**\
The referenced link may be anywhere (here at the bottom of the file)

![Figure 2]

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

## Links

**Absolute:** [Python installation](/p/Installation.md)\
**Relative:**[Matlab installation](../m/Installation.md)\
**Image** [fig1](/assets/images/output_69_0.png)

**link**: {%link p/Installation.md %}\
The link tag cannot be used on github pages because it still does not prepend
site.baseurl.

[Figure 2]: {{ "assets/images/output_69_0.png" | relative_url }}