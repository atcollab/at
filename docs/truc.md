---
title: Test
pict2: /assets/images/output_69_0.png
---
{% assign pict1 = "/assets/images/output_69_0.png" | relative_url %}
{% assign pict2 = "/assets/images/output_69_0.png" | absolute_url %}

base_url: {{ base_url }}

absolute_img: {{ page.pict2 | absolute_url}}

relative_img: {{ page.pict2 | relative_url}}

abs pict: {{ pict2 }}

rel pict: {{ pict1 }}

{%link assets/images/output_69_0.png %}
{%include image.html file="/assets/images/output_69_0.png" caption="Figure 1: test caption" %}
