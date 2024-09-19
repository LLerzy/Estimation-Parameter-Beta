# Estimación Posterior De Los Parámetros De Forma De La Distribución Beta: Estudio de Simulación

En este repositorio se presentan un conjunto de funciones diseñadas para estimar desde el enfoque Bayesiano, los parámetros de forma de la distribución Beta, cuya función de densidad es

$$f_X(x|\alpha,\beta)=\dfrac{1}{B(\alpha,\beta)}x^{\alpha-1}(1-x)^{\beta-1};\hspace{1cm}x\in[0,1],$$

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
    Rscript Post_Estimation_Beta_Param.R
    ```

El código considera tres escenarios para los parámetros de forma de la distribución beta y establece diferentes sesgos para el grado de información de expertos hipotéticos. Los resultados que genera el código, involucra el enfoque empírico de Bayes y el enfoque subjetivo, dentro de estos se presenta: Matriz de valores de hiperparámetros, matriz de estimaciones (teóricas) de las distribuciones a priori, gráficos sobre las estimaciones posteriores generadas por el método de muestreo por importancia y características de los estimadores posteriores (estimación promedio, sesgo, error cuadrático medio, probabilidad de cobertura, longitud promedio).

## Archivos y Documentos

-   [`Post_Estimation_Beta_Param.R`](Post_Estimation_Beta_Param.R): Contiene el código utilizado para obtener los resultados del estudio de simulación sobre la estimación posterior de los parámetros de forma de la distribución beta, considerando una nueva distribución bivariada como a priori y una propuesta metodológica para la obtención de valores de hiperparámetros.
-   [`requiredfunctions.R`](requiredfunctions.R): Contiene las funciones personalizadas utilizadas en el estudio de simulación.
-   [`Results_Post_Estim/`](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/Post-Estimate/Results_Post_Estim): Carpeta que contiene dos subcarpetas, la primera almacena las figuras generadas para el monitoreo de las características de los estimadores posteriores y la segunda contiene hojas de calculo sobre los valores de hiperparámetros, junto con las cadenas utilizadas en la construcción de los gráficos antes descritos.

### Figuras y Documentos Claves:

-   **Estimaciones Posteriores Generadas Para El Escenario $\alpha=0.5$ y $\beta=0.5$**: Consulte la figura generada [`Alpha-Escenario1.png`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Post-Estimate/Results_Post_Estim/Graphics/ParAlphaSigMu01SigV01Scen-1.png) para los resultados obtenidos del parámetros $\alpha$ y la figura [`Beta-Escenario1.png`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Post-Estimate/Results_Post_Estim/Graphics/ParBetaSigMu01SigV01Scen-1.png) para los resultados de $\beta$.
-   **Suplemento 1**: Es un archivo PDF, [`Max-IC-Varianza`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Post-Estimate/Results_Post_Estim/Supplement1.pdf), que contiene los resultados obtenidos para las estimaciones posteriores cuando se considera que la varianza se encuentra acotada por el límite superior del intervalo de cuantiles que proporciona el experto o el enfoque emppirico de Bayes.

## Contribuciones

¡Las contribuciones son bienvenidas! Envíe una solicitud de incorporación de cambios o abra un problema si tiene alguna sugerencia o mejora.
