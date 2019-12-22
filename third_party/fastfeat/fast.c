// clang-format off
#include <stdlib.h>
#include "fast.h"


xy* eb_aom_fast9_detect_nonmax(const byte* im, int xsize, int ysize, int stride, int b, int* ret_num_corners)
{
    xy* corners;
    int num_corners;
    int* scores;
    xy* nonmax;

    corners = eb_aom_fast9_detect(im, xsize, ysize, stride, b, &num_corners);
    scores = eb_aom_fast9_score(im, stride, corners, num_corners, b);
    nonmax = eb_aom_nonmax_suppression(corners, scores, num_corners, ret_num_corners);

    free(corners);
    free(scores);

    return nonmax;
}
// clang-format on
