![image](https://github.com/TheHopeSharedFoundation/ReLiEF-fingerprints/assets/7951822/e459a3ff-68ac-4e5b-b9a6-024f425b246b) 

# ReLiEF fingerprints

Receptor Ligand surfacE Feature (ReLiEF) fingerprints for encoding molecular structure and physicochemical properties.  Essentially, ReLiEF fingerprints are generated from a molecular surface that is colored by a particular property, such as partial charge.  These are physicochemical 3D-based fingerprints and can be used in a variety of applications such as similarity clustering, machine learning, and pharmacophore generation.  There are two types of ReLiEF fingerprints, segmentation-based and distribution-based.

### Segmentation-based ReLiEF fingerprints

In segmentation-based ReLiEF fingerprints (sb-ReLiEF), the molecule is first oriented with respect to a common coordinate system.  The whole surface (continous surface) of the molecule is captured in a VRML file with vertex points and associated colors.  Then the surface is sliced into equally-spaced segments.  A band is moved across each segment and the points of intersection of the segment with the band are recorded as fingerprint bits.  The bits can simultaneously capture a variety of information such as distance from a plane, pharmacophore, and atom type.  The bits from each segment are then concatenated into a fingerprint.  The main requirement of sb-ReLiEF is that all structures be aligned to a common coordinate system.  This video describes the use of sb-ReLiEF and TMAP (https://tmap.gdb.tools/) to cluster and visualize the structural similarity of Abl kinase conformations (from [DOI: 10.1126/science.abc2754](https://www.science.org/doi/10.1126/science.abc2754)): https://www.youtube.com/watch?v=2tnHeRlLt7Y 

### Distribution-based ReLiEF fingerprints

In contrast to sb-ReLiEF, distribution-based ReLiEF fingerprints (db-ReLiEF) do not depend on molecular orientation or position.  In db-ReLiEF, the surface of the molecule is first represented as a dot-surface (discrete surface) colored by a particular property and captured in a VRML file.  The dots representing the surface are then binned according to their color.  For each color bin, the property that the color represents and the sequence of dot addition to the bin (index/instance) are recorded with each dot, and this represents a fingerprint bit.  This binning process results in a distribution of the physicochemical 3D-based property represented by the color.  This process can be repeated for multiple properties, generating a distribution for each property.  For each property distribution and each bin, the fingerprint bits are concatenated in order of dot index into a fingerprint.  In this way, the fingerprint captures the distributions of multiple properties simultaneously, such as partial charge, distance from the center of mass, and functional group.  This schematic describes the process of creating db-ReLiEF based on the partial charges of PDE-10 ligands (from [DOI:10.1007/s10822-022-00478-x](https://pubmed.ncbi.nlm.nih.gov/36153472/)):        


![Screen Shot 2023-05-30 at 10 53 32 AM](https://github.com/TheHopeSharedFoundation/ReLiEF-fingerprints/assets/7951822/c54aebd0-6815-4325-a57e-edeb76e0daff)


![Screen Shot 2023-05-30 at 9 02 44 AM](https://github.com/TheHopeSharedFoundation/ReLiEF-fingerprints/assets/7951822/fe4200a3-d66d-41d3-9c5b-64083dae16b1)

This is a repository and collection of work by the Hope Shared Foundation, an all-volunteer 501(c)(3) not-for-profit dedicated to empowering students, citizen-scientists, and researchers to contribute to open-access biomedical research and innovation.  All work in this repository is licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
