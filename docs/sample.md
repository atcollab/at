---
title: Test
pict2: /assets/images/output_69_0.png
---
{% assign pict1 = "/assets/images/output_69_0.png" %}

## Heading with explicit id {#location}

## Heading with specified class
{:.red}

**Using variables:**

base_url: {{ site.base_url }}

absolute_img: {{ page.pict2 | absolute_url}}

relative_img: {{ page.pict2 | relative_url}}

abs pict: {{ pict1 | absolute_url }}

rel pict: {{ pict1 | relative_url }}

## Heading with auto id

{%link assets/images/output_33_1.png %}

**In-line formula:** $$\eta_c = 1/\gamma^2 - \alpha_c$$

**Display formula:**

$$\gamma=\frac{1+\alpha^2}{\beta}$$



**Image with in-line link:**

![Figure 1]({{ "/assets/images/output_33_1.png" | relative_url}})

**Image with referenced link:**\
The referenced link may be anywhere (here at the bottom of the file)

![Figure 2]

**Image with format and caption** (specific for the AT theme):
{%include image.html src="/assets/images/output_67_1.png" caption="Figure 3: test caption" %}

[Figure 2]: {{ "/assets/images/output_69_0.png" | relative_url}}