

// store the mesh index near IB points
std::vector<std::vector<label>> meshSupport;
// store the processor number near IB points
std::vector<std::vector<label>> meshSupportProcessor;
// store the Dirac delta function value near IB points
std::vector<std::vector<scalar>> Dirac;

//const Foam::fvMesh& B_mesh = U.mesh();

//mesh spacing
//const scalar h = 0.05;
scalar h = GREAT;
forAll(mesh.C(), cellI){
    const cell& c = mesh.cells()[cellI];
    scalar minEdgeLength = GREAT;
    forAll(c, faceI) {
        const face& f = mesh.faces()[c[faceI]];
        forAll(f, pointI) {
        minEdgeLength = min(minEdgeLength, mag(mesh.points()[f[pointI]] - mesh.points()[f.nextLabel(pointI)]));
        }
    }
    h = min(h, minEdgeLength);
}
Info << "Minimum cell length in the domain is: " << h << endl;

/*
const edgeList& IBedges = mesh.edges();
const pointField& IBPoint = mesh.points();

// Make a list with all the edge lenghts
scalarField MeshLengths(IBedges.size(), 0.0);

forAll(IBedges, EdgeI)
{
    MeshLengths[EdgeI] = IBedges[EdgeI].mag(IBPoint);
}

// Find the minimum edge length in the whole mesh
scalar minMeshSize = Foam::min(MeshLengths);
scalar h = minMeshSize;

Info << "The minimum edge length in the whole mesh: " << h << endl;
 */

// npt is the number of IB points
for(label I=0; I<npt; I++)
{
    scalar X = Lmarks[I].x();
    scalar Y = Lmarks[I].y();

    std::vector<label> Temp1;
    meshSupport.push_back(Temp1);
    meshSupportProcessor.push_back(Temp1);
    
    std::vector<scalar> Temp2;
    Dirac.push_back(Temp2);
    // loop over all the processor
	for(int ip=0;ip<mpi_num;ip++)
	{
        // loop over all the cells in the processor
	    for(int J=0;J<cellNum[ip];J++)
	    {
        //cal the distance between the IB points and the mesh centers

        scalar rx = fabs((meshCenterLocal_Collect[0][ip][J] - X)/h);
		scalar ry = fabs((meshCenterLocal_Collect[1][ip][J] - Y)/h);

        //use rxc and ryc determine whether it is within the scope of influence
		scalar rxc = rx - 1.5;
		scalar ryc = ry - 1.5;

		if ((rxc < 0.0001) && (ryc < 0.0001))
		{
            //store the mesh index near IB points
		    meshSupport[I].push_back(J);
            //store the processor number near IB points
            meshSupportProcessor[I].push_back(ip);
		    scalar Faix = Delta(rx);
		    scalar Faiy = Delta(ry);

		    Dirac[I].push_back(Faix*Faiy/h/h);      
		}
	    }
	}
}
if(myid==0)
{
    Info << "Number of nearby cells = " << meshSupport[0].size() << endl;
}
MPI_Barrier(MPI_COMM_WORLD);
//Info << "Lmarks[0] coordinates: " << "X = " << Lmarks[0].x() << "  " << "Y = " << Lmarks[0].y() << endl;


/*
label nptZero = meshSupport[0].size();
for(label I=0; I<nptZero; I++)
{
    label NN = meshSupport[0][I];
    label NP = meshSupportProcessor[0][I];
    Info << "No." << I+1 << ": " << NN << endl;
    Info << "Coordinates: " << "x = " << meshCenterLocal_Collect[0][NP][NN] << "  " << "y = " << meshCenterLocal_Collect[1][NP][NN] << "  " << "z = " << meshCenterLocal_Collect[2][NP][NN] << endl;
    Info << "x-X: " << meshCenterLocal_Collect[0][NP][NN] - Lmarks[0].x() << "   " << "y-Y: " << meshCenterLocal_Collect[1][NP][NN] - Lmarks[0].y() << endl;
    Info << "Dirac = " << Dirac[0][I] << endl;
}
*/

//Info << "Nearest Cell - 136481: x:" << meshCenters[136481].x() << "   " << "y: " << meshCenters[136481].y() << endl;

/*
for(int I=0; I<m; I++)
{
    int n = meshSupport[I].size();
    for(int J=0; J<n; J++)
    {
        Info<< meshSupport[I][J] << "   ";
    }
    printf("\n");
}
*/
    
