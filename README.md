# Estimación Posterior De Los Parámetros De Forma De La Distribución Beta: Estudio de Simulación

Este repositorio es una rama de uno principal denominado [**Estimation-Parameter-Beta**](https://github.com/LLerzy/Estimation-Parameter-Beta).

## Tabla de Contenido

-   [Resumen](#resumen)
-   [Instalación](#instalación)
-   [Uso](#uso)
-   [Archivos y Documentos](#archivos-y-documentos)
-   [Contribuciones](#contribuciones)

## Resumen {#resumen}

Este proyecto contiene un conjunto de rutinas destinadas a estimar los parámetros de forma de la distribución Beta para la variable $X$ desde una perspectiva bayesiana. Aunque se utiliza una nueva distribución a priori bivariada y su respectivo método para simular muestras aleatorias (los cuales pueden ser consultados en la rama ), el objetivo principal es desarrollar un estudio de simulación para monitorear el comportamiento de las estimaciones obtenidas para estos parámetros.

El script incluyen enfoques bayesianos empíricos y subjetivos para la estimación de hiperparámetros, utilizando tanto intervalos bootstrap como intervalos de expertos establecidos con diferentes sesgos y semi-amplitudes.

Todo el código está escrito en R-Project, aprovechando las herramientas estadísticas y computacionales de este software. El proyecto incluye tanto resultados numéricos como teóricos.

## Instalación {#instalación}

Para ejecutar el código es necesario tener `R` instalado con los siguiente paquetes: - `ggplot2` - `gridExtra` - `tidyr` - `plotly` - `coda` - `foreach` - `doParallel` - `betafunctions` - `openxlsx` - `xtable`

Usted puede instalar los paquetes requeridos utilizando el siguiente comando:

``` r
install.packages(c("ggplot2", "gridExtra", "tidyr","plotly","coda","foreach","doParallel","betafunctions","openxlsx","xtable"))
```

## Uso {#uso}

Para repliar el análisis y ejecutar los algoritmos, siga los siguientes pasos:

1.  Clone este repositorio:

    ``` bash
    git clone https://github.com/LLerzy/Estimation-Parameter-Beta/tree/Algorithm-Sim-Samples.git
    ```

2.  Abra el script de R o archivo RMarkdown en RStudio o en cualquier ambiente de desarrollo integrado (IDE) compatible:

    -   Para el script que contiene las funciones principales, consulte [`requiredfunctions.R`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Algorithm-Sim-Samples/requiredfunctions.R).
    -   Alternativamente, puede ver la salida de Markdown renderizada en [`randomsamples.md`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Algorithm-Sim-Samples/randomsamples.md).

3.  Corra el script para generar los resultados:

    ``` bash
    Rscript -e "rmarkdown::render('randomsamples.Rmd')"
    ```

El código generará los diagnósticos de convergencia, como histograma versus densidad, gráficos de calor para el tamaño de muestra efectivo, tasa de aceptaicón y R-Hat, promedio acumulado, traza, autocorrelación y curvas de nivel junto con diagrama de puntos de la muestra generada.

## Archivos y Documentos {#archivos-y-documentos}

-   [`randomsamples.Rmd`](randomsamples.Rmd): El archivo RMarkdown principal que contiene el código para generar los resultados.
-   [`randomsamples.md`](randomsamples.md): Una versión Markdown de los resultados.
-   [`requiredfunctions.R`](requiredfunctions.R): Contiene las funciones personalizadas utilizadas en el algoritmo, como el Gibbs Sampling, Metropolis Hasting y la nueva distribución bivariada.
-   [`randomsamples_files/`](randomsamples_files/figure-gfm): Carpeta que contiene las figuras generadas por el archivo RMarkdown.

### Figuras Clave:

-   **Tamaño Efectivo y Tasa de Aceptación según Precisión vs. Tamaño de Muestra**: Consulte la figura generada [`EffectiveSize-AcceptanceRate.png`](randomsamples_files/figure-gfm/unnamed-chunk-3-1.png) para analizar el rendimiento del algoritmo en diferentes configuraciones de la precisión y el tamaño de muestra.
-   **Curvas de Nivel vs Diagrama de Puntos del Vector Aleatorio** $(Y_1,Y_2)$: Consulte el gráfico en [`Curvas-de-Nivel-vs-Diagrama-Puntos.png`](randomsamples_files/figure-gfm/unnamed-chunk-14-1.png) para obtener una evaluación de la convergencia de las cadenas.

## Contribuciones {#contribuciones}

¡Las contribuciones son bienvenidas! Envíe una solicitud de incorporación de cambios o abra un problema si tiene alguna sugerencia o mejora.
