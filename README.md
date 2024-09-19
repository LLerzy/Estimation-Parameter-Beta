# Estimación Posterior De Los Parámetros De Forma De La Distribución Beta: Estudio de Simulación

Este repositorio principal contiene tres ramas, las cuales tienen como objetivo proporcionar una metodología para estimar los parámetros de forma de la distribución Beta.

## Tabla de Contenido

-   [Resumen](#resumen)
-   [Instalación](#instalación)
-   [Uso](#uso)
-   [Funciones Incluidas](#funciones-incluidas)
-   [Ramas](#ramas)
-   [Contribuciones](#contribuciones)

## Resumen

En este repositorio se presentan un conjunto de funciones diseñadas para estimar desde el enfoque Bayesiano, los parámetros de forma de la distribución Beta, cuya función de densidad es

$$f_X(x|\alpha,\beta)=\dfrac{1}{beta(\alpha,\beta)}x^{\alpha-1}(1-x)^{\beta-1};\hspace{1cm}x\in[0,1],$$

siendo $beta(\alpha,\beta)$ la función beta y $\alpha>0$ y $\beta>0$. Para los parámetros de forma $(\alpha,\beta)$ se utiliza una nueva distribución bivariada para la cual se han encontrado de forma analítica sus momentos conjuntos, se ha diseñado un algoritmo de simulación de muestras aleatorias, se ha monitoreado la convergencia de diferentes cadenas generadas por medio de dicho algoritmo, y se ha desarrollado un estudio de simulación considerando diferentes escenarios para los parámetros de forma antes mencionados.

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

3.  Puede ejecutar el script y cargar todas las funciónes con el siguiente comando: 

    ``` bash
    Rscript requiredfunctions.R
    ```

## Funciones Incluidas

Las funciones que son definidas en el script `requiredfunctions.R` son las siguientes:

-   `Prior`: Define la función de densidad de probabilidad a priori propuesta denominada combinación 1 en la rama [New-Biv-Dist](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/New-Biv-Dist).
-   `FC_X1_Given_v`: Distribución condicional completa de $X_1$ dado $X_2 = v$.
-   `Graph_Fc_X1`: Grafica la distribución condicional completa para los valores de parámetros dados.
-   `Gen_FC_X1_X2`: Algoritmo Metropolis-Hastings con caminata aleatoria para generar muestras de la distribución condicional completa.
-   `Mon_Measure`: Monitorea las tasas de aceptación y el ESS (tamaño efectivo de la muestra) para diferentes valores de $v$ y precisión.
-   `Mon_R_Hat`: Monitorea el diagnóstico Gelman-Rubin ($R$-hat) para diferentes valores de $v$ y precisión.
-   `Graphs`: Genera gráficos de histogramas, densidad, trazas y control de convergencia utilizando el promedio.
-   `Gen_Joint_Dist`: Gibbs Sampling para generar muestras de la distribución conjunta del vector aleatorio $(X_1,X_2)$.
-   `Mtovar_vs2`: Generaliza el método de Tovar para obtener los valores de los hiperparámetros.
-   `Mom_Prior_Dist`: Calcula los momentos conjuntos de orden $l$ para la distribución a priori propuesta.
-   `Measure_Diagnostic`: Compara los resultados analíticos y numéricos para muestras de datos proporcionados por el usuario.
-   `Measure_Analy`: Calcula los resultados analíticos para la distribución a priori propuesta.
-   `Hyperparameters`: Obtiene los valores de los hiperparámetros utilizando enfoques Bayesianos empírico y subjetivo.
-   `Est_Post`: Estimación posterior de los parámetros $\alpha$ y $\beta$ de la distribución Beta utilizando muestreo por importancia.
-   `Sim_study`: Realiza estudios de simulación para comparar estimaciones posteriores utilizando diferentes hiperparámetros y tamaños de muestra.
-   `Individual_Graphs`: Crea gráficos individuales para monitorear las estimaciones posteriores utilizando el sesgo y el MSE como indicadores.
-   `Comparison_Hyper`: Compara las funciones de densidad conjuntas para diferentes conjuntos de hiperparámetros.

### Notas:

-   Se recomienda revisar y adaptar cada función según las necesidades específicas de cada análisis.
-   Asegúrate de comprender cada función antes de usarla para garantizar resultados precisos y evitar errores potenciales.

## Ramas

En este repositorio se encuentran las siguientes ramas:

-   [**New-Biv-Dist**](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/New-Biv-Dist): Presenta características aproximadas numéricamente, para siete de las 25 distribuciones de probabilidad propuestas para representar el comportamiento aleatorio de los parámetros de forma de la distribución Beta. Además, se comparan algunas aproximaciones numéricas y teóricas de la distribución denominada combinación 1.
-   [**Algorithm-Sim-Samples**](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/Algorithm-Sim-Samples): Monitoreo de la convergencia de las cadenas generadas por un algoritmo que utiliza métodos MCMC para generar muestras aleatorias de una nueva distribución de probabilidad propuesta, denominada combinación 1 en la rama [New-Biv-Dist](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/New-Biv-Dist).
-   [**Post-Estimate**](https://github.com/LLerzy/Estimation-Parameter-Beta/tree/Post-Estimate): Desarrolla un estudio de simulación para estimar los parámetros de forma de la distribución Beta. Este estudio involucra: dos preguntas que pueden ser utilizadas en un proceso de elicitación, método de obtención de hiperparámetros, tres escenarios para los parámetros de forma, información de expertos hipotéticos con diferentes niveles de sesgo y gráficas que comparan estimadores posteriores a través de la estimación promedio, sesgo, error cuadrático medio, probabildiad de cobertura y longitud promedio.

## Contribuciones

¡Las contribuciones son bienvenidas! Envíe una solicitud de incorporación de cambios o abra un problema si tiene alguna sugerencia o mejora.
