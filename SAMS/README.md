# SAMS
SAMS: Surface Analysis, Mapping, and Statistics

This is a package containing software that allows for simultaneous matching of surfaces as well as a degree of statistical analysis. This package is suitable for rigid, homeomorphic, simply connected surfaces like anatomical surfaces and not for other datasets of interests such as SCAPE, FAUST, or non-homeomorphic or non-simply connected data.

***Instructions***

1. Download SAMS into your desired directory
2. Edit the following files as needed: ./Setup/AlignmentSetup.m, ./Setup/MappingSetup.m, ./Setup/StatisticsSetup.m
3. OPTIONAL: If data is not already aligned, run Execute00
4. Run the files starting in  Execute01, Execute02, and Execute03.

***Citations***

**If using this software for mapping or statistics, please cite the following works:**

Ravier, Robert J. *Algorithms with Applications to Anthropology.* Diss. Duke University, 2018.
Ravier, Robert J. "Eyes on the Prize: Improved Registration via Forward Propagation." *arXiv preprint arXiv:1812.10592 (2018).*
Ravier, Robert J et al. "SAMS: A New Method for Surface Alignment, Matching, and Statistics." *In preparation, 2019.*

**If you use the mapping package for topological spheres, please cite:**

Aigerman, Noam, and Yaron Lipman. "Hyperbolic orbifold tutte embeddings." *ACM Transactions on Graphics (TOG)* 35.6 (2016): 217.

The alignment package Auto3dgm (http://github.com/trgao10/PuenteAlignment) is included in the Alignment folder for convenience, though was modified so as to work with file structure and does not have the capabilities for distributed computing at the original package. **If you utilize Auto3dgm, please cite:**

Boyer, Doug M., et al. "A New Fully Automated Approach for Aligning and Comparing Shapes." *The Anatomical Record* 298.1 (2015): 249-276.

Puente, Jesús. *Distances and Algorithms to Compare Sets of Shapes for Automated Biological Morphometrics.* PhD Thesis, Princeton University, 2013.