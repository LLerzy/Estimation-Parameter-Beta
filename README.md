# Estimación Posterior De Los Parámetros De Forma De La Distribución Beta: Estudio de Simulación

Este repositorio principal contiene tres ramas, las cuales tienen como objetivo proporcionar una metodología para estimar los parámetros de forma de la distribución Beta.

## Tabla de Contenido

-   [Resumen](#resumen)
-   [Instalación](#instalación)
-   [Uso](#uso)
-   [Funciones Incluidas](#funciones-incluidas)
-   [Contribuciones](#contribuciones)

## Resumen

En este repositorio se presentan un conjunto de funciones diseñadas para estimar desde el enfoque Bayesiano, los parámetros de forma de la distribución Beta, cuya función de densidad es

$$f_X(x|\alpha,\beta)=\dfrac{1}{beta(\alpha,\beta)}x^{\alpha-1}(1-x)^{\beta-1};\hspace{1cm}x\in[0,1],$$

siendo $beta(\alpha,\beta)$ la función beta y $\alpha>0$ y $\beta>0$. Para los parámetros de forma $(\alpha,\beta)$ se utiliza una nueva distribución bivariada para la cual se han encontrado de forma analítica sus momentos conjuntos, se ha diseñado un algoritmo de simulación de muestras aleatorias, se ha monitoreado la convergencia de diferentes cadenas generadas por medio de dicho algoritmo, y se ha desarrollado un estudio de simulación considerando diferentes escenarios para los valores de parámetros de forma antes mencionados.

## Instalación

Para ejecutar el código es necesario tener `R` instalado con los siguiente paquetes: 

- `ggplot2` 
- `gridExtra` 
- `tidyr` 
- `plotly` 
- `coda` 
- `foreach` 
- `doParallel` 
- `betafunctions` 
- `openxlsx` 
- `xtable`

Usted puede instalar los paquetes requeridos utilizando el siguiente comando:

``` r
install.packages(c("ggplot2", "gridExtra", "tidyr","plotly","coda","foreach","doParallel","betafunctions","openxlsx","xtable"))
```

## Uso

Para poder ejecutar los algoritmos diseñados, siga los siguientes pasos:

1.  Clone este repositorio:

    ``` bash
    git clone https://github.com/LLerzy/Estimation-Parameter-Beta/tree/main.git
    ```

2.  Abra el script de R en cualquier ambiente de desarrollo integrado (IDE) compatible:

    -   Para el script que contiene las funciones principales en el proceso de estimación de los parámetros de forma, consulte [`requiredfunctions.R`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Post-Estimate/requiredfunctions.R).

3.  Puede ejecutar el script y carargar todas las funciónes con el siguiente comando, 

    ``` bash
    Rscript requiredfunctions.R
    ```

El código considera tres escenarios para los parámetros de forma de la distribución beta y establece diferentes sesgos para el grado de información de expertos hipotéticos. Los resultados que genera el código, involucra el enfoque empírico de Bayes y el enfoque subjetivo, dentro de estos se presenta: Matriz de valores de hiperparámetros, matriz de estimaciones (teóricas) de las distribuciones a priori, gráficos sobre las estimaciones posteriores generadas por el método de muestreo por importancia y características de los estimadores posteriores (estimación promedio, sesgo, error cuadrático medio, probabilidad de cobertura, longitud promedio).

## Funciones Incluidas

-   [`Post_Estimation_Beta_Param.R`](Post_Estimation_Beta_Param.R): Contiene el código utilizado para obtener los resultados del estudio de simulación sobre la estimación posterior de los parámetros de forma de la distribución beta, considerando una nueva distribución bivariada como a priori y una propuesta metodológica para la obtención de valores de hiperparámetros.

-   `Prior`: Define la función de densidad de probabilidad a priori propuesta.
-   `FC_X1_Given_v`: Distribución condicional completa de \(X_1\) dado \(X_2 = v\).
-   `Graph_Fc_X1`: Grafica la distribución condicional completa para los valores de parámetros dados.
-   `Gen_FC_X1_X2`: Algoritmo Metropolis-Hastings con caminatas aleatorias para generar muestras de la distribución condicional completa.
-   `Mon_Measure`: Monitorea las tasas de aceptación y el ESS (tamaño efectivo de la muestra) para diferentes valores de \(v\) y precisión.
-   `Mon_R_Hat`: Monitorea el diagnóstico Gelman-Rubin (\(R\)-hat) para diferentes valores de \(v\) y precisión.
-   `Graphs`: Genera gráficos de histogramas, densidad, trazas y control de convergencia utilizando el promedio.
-   `Gen_Joint_Dist`: Muestreo Gibbs para generar distribuciones conjuntas de \(X_1\) y \(X_2\).
-   `Mtovar_vs2`: Generaliza el método de Tovar para obtener los valores de los hiperparámetros.
-   `Mom_Prior_Dist`: Calcula los momentos conjuntos de orden \(l\) para la distribución a priori propuesta.
-   `Measure_Diagnostic`: Compara los resultados analíticos y numéricos para muestras de datos dadas.
-   `Measure_Analy`: Calcula los resultados analíticos para la distribución a priori propuesta.
-   `Hyperparameters`: Obtiene los hiperparámetros utilizando enfoques Bayesianos empírico y subjetivo.
-   `Est_Post`: Estimación posterior de los parámetros \(\alpha\) y \(\beta\) de la distribución Beta utilizando muestreo por importancia.
-   `Sim_study`: Realiza estudios de simulación para comparar estimaciones posteriores utilizando diferentes hiperparámetros y tamaños de muestra.
-   `Individual_Graphs`: Crea gráficos individuales para monitorear las estimaciones posteriores utilizando el sesgo y el MSE como indicadores.
-   `Comparison_Hyper`: Compara las funciones conjuntas para diferentes conjuntos de hiperparámetros.

### Notas:

-   Se recomienda revisar y adaptar cada función según las necesidades específicas de cada análisis.
-   Asegúrate de comprender cada función antes de usarla para garantizar resultados precisos y evitar errores potenciales.

## Contribuciones

¡Las contribuciones son bienvenidas! Envíe una solicitud de incorporación de cambios o abra un problema si tiene alguna sugerencia o mejora.
