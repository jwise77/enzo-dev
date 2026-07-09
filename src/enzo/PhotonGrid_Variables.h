#define MAX_SOURCES 3000
// Photons
int    NumberOfPhotonPackages;
//int    NumberOfRenderingPackages;

//  SoA of active packages
PhotonPackageSoA PhotonPackages;

// SoA of packages with its work already finished
PhotonPackageSoA FinishedPhotonPackages;  

// SoA of packages that are "paused", waiting to be merged
PhotonPackageSoA PausedPhotonPackages;


// linked list of packages used for projections or volume renderings
//PhotonPackageEntry *RenderingPhotonPackages;

// Sources
//int    NumberOfRadiationSources;
grid  **SubgridMarker; // pointers to first array elements of subgrids for each cell
int HasRadiation;

// For adaptive timestep control while restricting the change in HII
// to 50%, we need to record the minimum kph for which a ray passes
// through with a cumulative optical depth >0.1.
float MaximumkphIfront;
int IndexOfMaximumkph;
int OriginalProcessorNumber;
