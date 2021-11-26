---
title: Test
pict2: /assets/images/output_69_0.png
---
{% assign pict1 = "/assets/images/output_69_0.png" %}
var1: {{ pict1}}

var2: {{ page.pict2}}

{% link assets/images/output_69_0.png %}
{%include image.html file="/assets/images/output_69_0.png" caption="Figure 1: test caption" %}
