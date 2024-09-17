# Algoritmo Basado En Métodos MCMC Para Simular Muestras De Una Nueva Distribución Bivariada

Este repositorio es una rama del repositorio principal denominado Estimation-Parameter-Beta.

## Tabla de Contenido
- [Resumen](#resumen)
- [Instalación](#instalación)
- [Uso](#uso)
- [Archivos y Documentos](#archivos-y-documentos)
- [Contributing](#contributing)
- [License](#license)

## Resumen

Este proyecto presenta un conjunto de funciones que implementan algoritmos MCMC diseñados para generar muestras aleatorias del vector bivariado \((X_1, X_2)\) y posteriomente transformar dichas muestras al espacio del vector aleatorio $(Y_1,Y_2)$. Todo el código está escrito en R, aprovechando las herramientas estadísticas y computacionales de este software. El proyecto incluye tanto resultados numéricos como teóricos, con un enfoque en el análisis de las cadenas generadas a través de ocho criterios diferentes de convergencia. Estos diagnósticos proporcionan una evaluación integral del rendimiento de las cadenas, asegurando la fiabilidad y precisión de las estimaciones. Los métodos presentados son adecuados para una amplia gama de aplicaciones que involucran inferencia bayesiana y técnicas de muestreo de Gibbs.

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
```R
install.packages(c("ggplot2", "gridExtra", "tidyr","plotly","coda","foreach","doParallel","betafunctions","openxlsx","xtable"))
```

## Uso
Para repliar el análisis y ejecutar los algoritmos, siga los siguientes pasos:

1. Clone este repositorio:
   ```bash
   git clone https://github.com/LLerzy/Estimation-Parameter-Beta/tree/Algorithm-Sim-Samples.git
   ```
2. Abra el script de R o archivo RMarkdown en RStudio o en cualquier ambiente de desarrollo integrado (IDE) compatible:
   - Para el script que contiene las funciones principales, consulte [`requiredfunctions.R`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Algorithm-Sim-Samples/requiredfunctions.R).
   - Alternativamente, puede ver la salida de Markdown renderizada en [`randomsamples.md`](https://github.com/LLerzy/Estimation-Parameter-Beta/blob/Algorithm-Sim-Samples/randomsamples.md).

3. Corra el script para generar los resultados:
   ```bash
   Rscript randomsamples.md
   ```

El código generará los diagnósticos de convergencia, como Histograma versus Densidad, gráficos de calor para el tamaño de muestra efectivo, tasa de aceptaicón y R-Hat, promedio acumulado, traza, autocorrelación y curvas de nivel junto con diagrama de puntos de la muestra generada.

## Archivos y Documentos
- [`randomsamples.Rmd`](randomsamples.Rmd): El archivo RMarkdown principal que contiene el código para generar los resultados.
- [`randomsamples.md`](randomsamples.md): Una versión Markdown de los resultados.
- [`requiredfunctions.R`](requiredfunctions.R): Contiene las funciones personalizadas utilizadas en el algoritmo, como el Gibbs Sampling, Metropolis Hasting y la nueva distribución
- [`randomsamples_files/`](randomsamples_files/figure-gfm): Carpeta que contiene las figuras generadas por el archivo RMarkdown.

### Key Figures:
- **Effective Size vs Acceptance Rate**: Refer to the generated figure [`EffectiveSize-AcceptanceRate.png`](randomsamples_files/figure-gfm/EffectiveSize-AcceptanceRate.png) to analyze how the algorithm performs under different configurations.
- **R-Hat Convergence Plot**: See the R-Hat plot in [`R-hat.png`](randomsamples_files/figure-gfm/R-hat.png) for an evaluation of the convergence of the chains.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue if you have any suggestions or improvements.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
