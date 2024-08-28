# H&E-DAPI registration

For each sample,
1. Cut the DAPI via ```seperate_samples_LEVEL1.ipynb``` and save TIFF and anndata objects.
2. Cut H&E images via Fiji.
3. Do the registarion with ```registration_final_elastix.ipynb```

As for visualization, we need to read the anndata for each sample, then attach the corresponding registared H&E image to it as in ```plot_example.ipynb```.


# Phenocycler-DAPI Registration

1. Cut the DAPI via ```seperate_samples_LEVEL1_slide.ipynb```
2. Cut Phenocycler image ```Xenium Slide_Scan1.qptiff``` via Fiji.
3. Do the registarion with ```registration_final_elastix_phenocycler.ipynb```


**Cut Phenocycler image:**

1. Plugins > Bio-formats > Bio-formats importer
2. Check #2
3. Rotate 90deg right x2
4. Select sample and crop
5. Stack to image
6. Save each
