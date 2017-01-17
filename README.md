# distortion-model-precision

Goal: check the numerical precision of distortion models.
Principle:
- Generate pairs of distorted/undistorted points with a model from lensfun;
- Estimate the parameters of a given model from the pairs;
- Measure the RMSE and max error on independent point pairs.

The image coordinates are in [-1,1]. A correct precision is 1e-5, a great one 1e-16.

The model database of lensfun is in folder db. The script scripts/extract_distort_from_xml.sh permits to extract all distortion models from lensfun. You can inject one line in the program:

./precision_analysis -d "\<distortion model=poly3 focal=210 k1=0.00324 /\>" poly 4,9

Usage: ./precision_analysis [options] model order,order,...
Options:
        -r, --reverse Find inverse model
        -d, --distort=ARG XML line of lensfun db
 list of models:  radial  radial_center  radial_odd  radial_center_odd  division  division_center  division_even  division_center_even  FOV  poly
