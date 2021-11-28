---
title: Test
pict2: assets/images/output_69_0.png
---
{% assign pict1 = "/assets/images/output_69_0.png" %}

## Heading with explicit id {#location}

**Using variables:**

base_url: {{ site.baseurl }}

## Markdown features
{: #idfn .red}

**In-line formula:** $$\eta_c = 1/\gamma^2 - \alpha_c$$

**Display formula:**

$$\gamma=\frac{1+\alpha^2}{\beta}$$

**link:**
[Python installation]({%link p/Installation.md %})

**Image with in-line link:**

![Figure 1]({% link assets/images/output_33_1.png %})

**Image with referenced link:**\
The referenced link may be anywhere (here at the bottom of the file)

![Figure 2]

## Theme features

**Image with format and caption** (specific for the AT theme):
{%include image.html src="/assets/images/output_67_1.png"
alt="Figure 3" caption="Figure 3: test caption" %}

**Alerts:**
{% include note.html content="Alerts are displayed in colour boxes" %}

**Comments:** (do not appear)
{% comment %} The following lines are the referenced links
and do not appear in the output{% endcomment %}

[Figure 2]: {% link assets/images/output_69_0.png %}