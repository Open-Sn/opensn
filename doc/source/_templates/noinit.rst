{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

{%+ if objtype == 'class' -%}
.. autoclass:: {{ objname }}
   :members:
   :inherited-members:
{% endif %}

{%- if objtype == 'function' -%}
.. autofunction:: {{ objname }}
{% endif %}
