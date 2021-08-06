


//x-sweep
//loop over k
if (face == 0):
      istart = GridStartIndex[0];
      iend   = GridEndIndex[0];
      jstart = 0;
      jend   = GridDimension[1] - 1;
      kstart  
//y-sweep
//loop over i     
if (face == 1):
      jstart = GridStartIndex[1];
      jend   = GridEndIndex[1];
      kstart = 0;
      kend   = GridDimension[2] - 1;
//x-sweep
//loop over j      
if (face == 2):
      kstart = GridStartIndex[2];
      kend   = GridEndIndex[2];
      istart = 0;
      iend   = GridDimension[0] - 1;

//    |
//    |
//    | 
//    V

// outermost loop star and end indicies
int dim3Start[3] = {0,0,0};
int dim3End[3]   = {GridDimension[2], GridDimension[0], GridDimension[1]};
// mid loop star and end indicies
int dim2start[3] = {0,0,0};
int dim2end[3]   = {GridDimension[1] - 1, GridDimension[2] - 1, GridDimension[0] - 1};
// innermost loop star and end indicies
int dim1start[3] = {GridStartIndex[0], GridStartIndex[1], GridStartIndex[2]};
int dim1end[3]   = {GridEndIndex[0],   GridEndIndex[1], GridEndIndex[2]};



if (face == 0):
      indexA = (dim3 * GridDimension[1] + dim2) * GridDimension[0] + dim1;
      indexB = (dim3 * GridDimension[1] + dim2) * GridDimension[0] + dim1+1;
if (face == 1):
      indexA = (dim2 * GridDimension[1] + dim1) * GridDimension[0] + dim3;
      indexB = (dim2 * GridDimension[1] + dim1+1) * GridDimension[0] + dim3;
if (face == 2):
      indexA = ((dim1) * GridDimension[1] + dim3) * GridDimension[0] + dim2;
      indexB = ((dim1+1) * GridDimension[1] + dim3) * GridDimension[0] + dim2;


