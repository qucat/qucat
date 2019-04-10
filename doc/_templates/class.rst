{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

    {% block attributes %}
    
    {% if attributes %}
    .. rubric:: Attributes

    .. autosummary::
        :toctree: {{ objname }}
    {% for item in attributes %}
        ~{{ name }}.{{ item }}
    {%- endfor %}
    {% endif %}
    {% endblock %}

    {% block methods %}

    {% if methods %}
    .. rubric:: Methods

    .. autosummary::
        :toctree: {{ objname }}
    {% for item in methods %}
    {%- if item not in inherited_members %}
        ~{{ name }}.{{ item }}
    {%- endif %}
    {%- endfor %}
    {% endif %}
    {% endblock %}

