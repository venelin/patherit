# patherit
Pathogen Trait Heritability

Pathogen traits, such as the virulence of an infection, can vary significantly between patients. A major challenge is to measure the extent to which genetic differences between infecting strains explain the observed variation of the trait. This is quantified by the trait's broad-sense heritability, H^2. The patherit R-package allows to compare various estimators of pathogen trait heritability in data simulated using the `toyepidemic` R-package, and in real epidemiological data. 


The following methods are available:

1. Direct heritability estimator, R_{adj}^2, for data with full genetic and phenotype knowledge.
2. Donor recipient regression slope, b, for transmission couple data;
3. Several phylogenetic methods including ANOVA in phylogenetic pairs (PP) and closest phylogenetic pairs (CPP), the phylogenetic heritability estimated using either PMM or POUMM. 

For using the package, read the package [User Guide](./vignettes/UserGuide.Rmd).


# License

The toyepidemic R-package is licensed under version 3 of the [GNU General Public License](http://www.gnu.org/licenses/.

