# Nuevas Distribuciones De Probabilidad Bivariadas

Este repositorio es una rama de uno principal denominado [**Estimation-Parameter-Beta**](https://github.com/LLerzy/Estimation-Parameter-Beta).

## Tabla de Contenido

-   [Resumen](#resumen)
-   [Instalación](#instalación)
-   [Uso](#uso)
-   [Archivos y Documentos](#archivos-y-documentos)
-   [Contribuciones](#contribuciones)

## Resumen

Este proyecto presenta algunas distribuciones bivariadas construidas para el vector aleatorio $(Y_1,Y2)$ utilizando la transformación continua y diferenciable,  

$$\begin{matrix}
   T^{-1}: &  \mathbb{R}^2_+ & \Longrightarrow &(0,1)\times (0,U(X_1))\\
        &\left(Y_1,Y_2\right)&\longrightarrow & (X_1,X_2):=\left(\dfrac{Y_1}{Y_1+Y_2},\dfrac{Y_1Y_2}{(Y_1+Y_2)^2(Y_1+Y_2+1)}\right)
\end{matrix},$$

junto con la representación de la función de densidad conjunta del vector $(X_1,X_2)$,

$$f_{X_1,X_2}(x_1,x_2|\phi=(\phi_1,\phi_2)) = f_{X_1}(x_1|\phi_1)f_{X_2|X_1}(x_2|x_1,\phi_2).$$

Aunque se consideraron cinco distribuciones para $X_1$ y cinco para $X_2|X_1$, para un total de 25 combinaciones, esta rama únicamente presenta las caracteristicas de siete distribuciones, todas con vector de parámetros $\phi=(a,b,c,d)$. Para cada una, se construye la superficie de la densidad, curvas de nivel, densidad marginal, distribución de probabilidad acumulada marginal y momentos (media, varianza). Todo el código está escrito en Mathematica, aprovechando la manipulación simbólica y capacidad computacional de este software. El proyecto incluye resultados numéricos obtenidos con la función NIntegrate y algunos teóricos.

## Instalación

Para ejecutar el código es necesario tener `Mathematica` instalado.

## Uso

Para replicar el análisis y ejecutar los algoritmos, siga los siguientes pasos:

1.  Clone este repositorio:

    ``` bash
    git clone https://github.com/LLerzy/Estimation-Parameter-Beta/tree/New-Biv-Dist.git
    ```

2.  Abra el script de Mathematica en cualquier ambiente de desarrollo integrado (IDE) compatible:

    -   Para el script que contiene las características (mencionadas en el resumen) de las siete distribuciones bivariadas seleccionadas para el vector $(Y_1,Y_2)$, consulte [`SelectedBivariateDistributions.nb`](SelectedBivariateDistributions.nb).
    -   Para el script que contiene las caracteristicas de la distribución bivariada del vector $(Y_1,Y_2)$, construida considerando $X_1\sim Beta(a,b)$ y $X_2|X_1\sim Beta(0,U(X_1),c,d)$, consulte [`DistribucionBetaBeta4P.nb`](DistribucionBetaBeta4P.nb).

3.  El script debe ser ejecutado dentro del IDE.


El código en `DistribucionBetaBeta4P.nb` considera cuatro configuraciones para el vector de parámetros $(a,b,c,d)$ de la distribución bivariada construida utilizando $X_1\sim Beta(a,b)$ y $X_2|X_1\sim Beta(0,U(X_1),c,d)$, en el se presenta la superficie de la densidad, curvas de nivel, distribución de probabilidad acumulada marginal, superifice de supervivencia, función de densidad marginal y se comparan los momentos teóricos con los aproximados numéricamente a través de la función NIntegrate del software `Mathematica`.

## Archivos y Documentos

-   [`DistBivSelecc.pdf`](DistBivSelecc.pdf): Contiene resultados numéricos obtenidos para cada una de las siete distribuciones bivariadas seleccionadas.
-   [`DistribucionBetaBeta4P.pdf`](DistribucionBetaBeta4P.pdf): Contiene resultados numéricos y teóricos obtenidos para la distribución bivariada de $(Y_1,Y_2)$, construida al considerar $X_1\sim Beta(a,b)$ y $X_2|X_1\sim Beta(0,U(X_1),c,d)$.

## Contribuciones

¡Las contribuciones son bienvenidas! Envíe una solicitud de incorporación de cambios o abra un problema si tiene alguna sugerencia o mejora.
