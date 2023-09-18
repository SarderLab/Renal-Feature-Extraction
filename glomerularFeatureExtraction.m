function [features]=glomerularFeatureExtraction(segmentation_dir,image_dir)
% directory of glomerular component maps
segmented_gloms=dir([segmentation_dir,'/*.png']);
'Glomerular feature extraction'
Total=length(segmented_gloms);
features=zeros(Total,315);
% Remove small things smaller than this
min_object_size=50;

parfor q=1:Total
    features(q,:)=feature_extraction_inner_glomerulus(q,segmented_gloms,image_dir,min_object_size);
end