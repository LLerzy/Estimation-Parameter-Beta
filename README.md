# Estimación Posterior De Los Parámetros De Forma De La Distribución Beta: Estudio de Simulación

Este repositorio es una rama de uno principal denominado [**Estimation-Parameter-Beta**](https://github.com/LLerzy/Estimation-Parameter-Beta).

## Tabla de Contenido

-   [Resumen](#resumen)
-   [Instalación](#instalación)
-   [Uso](#uso)
-   [Archivos y Documentos](#archivos-y-documentos)
-   [Contribuciones](#contribuciones)

## Resumen

Este proyecto contiene un conjunto de rutinas destinadas a estimar los parámetros de forma de la distribución Beta para la variable $X$ desde una perspectiva bayesiana. Aunque se utiliza una nueva distribución a priori bivariada y su respectivo método para simular muestras aleatorias (los cuales pueden ser consultados en la rama [**Algorithm-Sim-Samples**](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/Algorithm-Sim-Samples)), el objetivo principal es desarrollar un estudio de simulación para monitorear el comportamiento de las estimaciones obtenidas para estos parámetros.

El script incluyen enfoques bayesianos empíricos y subjetivos para la estimación de hiperparámetros, utilizando tanto intervalos bootstrap como intervalos de expertos establecidos con diferentes sesgos y semi-amplitudes. Todo el código está escrito en R-Project, aprovechando las herramientas estadísticas y computacionales de este software. El proyecto incluye tanto resultados numéricos como teóricos.

## Instalación

Para ejecutar el código es necesario tener `R` instalado con los siguiente paquetes: - `ggplot2` - `gridExtra` - `tidyr` - `plotly` - `coda` - `foreach` - `doParallel` - `betafunctions` - `openxlsx` - `xtable`

Usted puede instalar los paquetes requeridos utilizando el siguiente comando:

``` r
install.packages(c("ggplot2", "gridExtra", "tidyr","plotly","coda","foreach","doParallel","betafunctions","openxlsx","xtable"))
```

## Uso

Para repliar el análisis y ejecutar los algoritmos, siga los siguientes pasos:

1.  Clone este repositorio:

    ``` bash
    git clone https://github.com/LLerzy/Estimation-Parameter-Beta/tree/Post-Estimate.git
    ```

2.  Abra el script de R en cualquier ambiente de desarrollo integrado (IDE) compatible:

    -   Para el script que contiene las funciones principales, consulte [`requiredfunctions.R`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Post-Estimate/requiredfunctions.R).
    -   Para el script que contiene las secuencias sobre el estudio de simulación, consulte [`Post_Estimation_Beta_Param.R`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Post-Estimate/Post_Estimation_Beta_Param.R).

3.  El script puede ser ejecutado con el siguiente comando, sin embargo, es necesario antes ajustar ciertos parámetros dentro del código.

    ``` bash
    Rscript -e "rmarkdown::render('randomsamples.Rmd')"
    ```

El código considera tres escenarios para los parámetros de forma de la distribución beta y establece diferentes sesgos para el grado de información de expertos hipotéticos. Los resultados que genera el código, involucra el enfoque empírico de Bayes y el enfoque subjetivo, dentro de estos se presenta: Matriz de valores de hiperparámetros, matriz de estimaciones (teóricas) de las distribuciones a priori, gráficos sobre las estimaciones posteriores generadas por el método de muestreo por importancia y características de los estimadores posteriores (estimación promedio, sesgo, error cuadrático medio, probabilidad de cobertura, longitud promedio).

## Archivos y Documentos

-   [`Post_Estimation_Beta_Param.R`](Post_Estimation_Beta_Param.R): Contiene el código utilizado para obtener los resultados del estudio de simulación para la estimación posterior de los parámetros de forma de la distribución beta, considerando una nueva distribución bivariada como a priori y una metodología de obtenciónd e valores de hiperparámetros.
-   [`requiredfunctions.R`](requiredfunctions.R): Contiene las funciones personalizadas utilizadas en el estudio de simulación.
-   [`Results_Post_Estim/`](randomsamples_files/figure-gfm): Carpeta que contiene las figuras generadas por el archivo RMarkdown.

### Figuras Clave:

-   **Tamaño Efectivo y Tasa de Aceptación según Precisión vs. Tamaño de Muestra**: Consulte la figura generada [`EffectiveSize-AcceptanceRate.png`](randomsamples_files/figure-gfm/unnamed-chunk-3-1.png) para analizar el rendimiento del algoritmo en diferentes configuraciones de la precisión y el tamaño de muestra.
-   **Curvas de Nivel vs Diagrama de Puntos del Vector Aleatorio** $(Y_1,Y_2)$: Consulte el gráfico en [`Curvas-de-Nivel-vs-Diagrama-Puntos.png`](randomsamples_files/figure-gfm/unnamed-chunk-14-1.png) para obtener una evaluación de la convergencia de las cadenas.

## Contribuciones

¡Las contribuciones son bienvenidas! Envíe una solicitud de incorporación de cambios o abra un problema si tiene alguna sugerencia o mejora.
