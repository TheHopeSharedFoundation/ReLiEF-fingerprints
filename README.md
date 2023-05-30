# ReLiEF fingerprints
Receptor Ligand surfacE Feature (ReLiEF) fingerprints for encoding molecular structure and physicochemical properties.  Essentially, ReLiEF fingerprints are generated from a molecular surface that is colored by a particular property, such as partial charge.  From the surface and color, fingerprints are generated.  These are physicochemical-3D-based fingerprints and can be used in a variety of applications such as similarity clustering, machine learning, and pharmacophore generation.  

There are two types of ReLiEF fingerprints, segmentation-based and distribution-based.  In segmentation-based ReLiEF fingerprints (sb-ReLiEF), the surface is first oriented with respect to a common coordinate system.  Then the surface is sliced into equally-spaced segments.  A band is moved across each segment and the points of intersection of the segment with the band are recorded as fingerprint bits.  The bits can simultaneously capture a variety of information such as distance from a plane, pharmacophore, and atom type.  The bits from each segment are then concatenated into a fingerprint.  The main advantage/disadvantage (depending on the particular application) of sb-ReLiEF is the requirement that all structures must be aligned to a common coordinate system.    




This video describes the use of ReLiEF fingerprints and TMAP (https://tmap.gdb.tools/) to cluster and visualize the structural similarity of Abl kinase conformations: https://www.youtube.com/watch?v=2tnHeRlLt7Y  


Click on the following Binder button to explore a TMAP tree-based heatmap visualization of ReLiEF fingerprints generated for Abl protein kinase: 

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/TheHopeSharedFoundation/ReLiEF-fingerprints/HEAD)

This is a repository and collection of work by the Hope Shared Foundation, an all-volunteer 501(c)(3) not-for-profit dedicated to empowering students, citizen-scientists, and researchers to contribute to open-access biomedical research and innovation.  All work in this repository is licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
