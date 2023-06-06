
Desde el Script "Ammonia_Water_Properties" importar la función PropsSI_NH3_H2O

Función que entrega las propiedades de la mezcla de agua con amoniaco.

Parámetros:
    - variable_required: Variable que será entregada por la función como resultado.
    - variable_1_name: Variable entregada como primer dato para obtener un resultado.
    - variable_1_value: Valor de la primera variable entregada como dato.
    - variable_2_name: Variable entregada como segundo dato para obtener un resultado.
    - variable_2_value: Valor de la segunda variable entregada como dato.

Retorna:
    - Valor de la variable solicitada en el primer parámetro, como resultado de las dos variables entregadas como datos.

Los parámetros 'variable_required', 'variable_1_name', y 'variable_2_name' pueden ser los siguientes:
    - 'p': Presión [Pa]
    - 'T': Temperatura [K]
    - 'Ql': Fracción másica de amoniaco en la fase líquida (entre 0 y 1)
    - 'Qg': Fracción másica de amoniaco en la fase gaseosa (entre 0 y 1)
    - 'Hl': Entalpía de la fase líquida [J/kg]
    - 'Hg': Entalpía de la fase gaseosa [J/kg]

    Los únicos pares de variables que no se pueden usar como variables de entrada son los siguientes:
        - 'p', 'Hl'
        - 'T', 'Hl'
        - 'Hl', 'Hg'
