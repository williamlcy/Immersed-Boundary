
//cal the distance between the first and second IB points
//delta s
scalar DS = Foam::sqrt((Lmarks[0].x()-Lmarks[1].x())*(Lmarks[0].x()-Lmarks[1].x()) + (Lmarks[0].y()-Lmarks[1].y())*(Lmarks[0].y()-Lmarks[1].y()));
//scalar DS = 0.0087;
Info << "Space of the Lmarks: " << DS << endl;
// face area between two IB points
scalar D_S = h*DS;

//the force of the structure to fluid
std::vector<std::vector<double>> F_Fluid_Collect[2];
F_Fluid_Collect[0].resize(mpi_num);
F_Fluid_Collect[1].resize(mpi_num);
for(int k=0;k<mpi_num;k++)
{
    F_Fluid_Collect[0][k].resize(cellNum[k],0.0);
    F_Fluid_Collect[1][k].resize(cellNum[k],0.0);
}

MPI_Barrier(MPI_COMM_WORLD);

//npt is the number of IB points
for(label I=0; I<npt; I++)
    {
        // the number of cells near the IB points
        label Nsupport = meshSupport[I].size();

        for(label J=0; J<Nsupport; J++)
        {
            // mesh cell global index
            label NNN = meshSupport[I][J];
            // mesh cell processor index
            label PPP = meshSupportProcessor[I][J];

            F_Fluid_Collect[0][PPP][NNN] += F_Structure[I].x()*Dirac[I][J]*D_S;
            F_Fluid_Collect[1][PPP][NNN] += F_Structure[I].y()*Dirac[I][J]*D_S;
        }
    }

//collect the force of the structure to fluid
forAll(F_Fluid, CELLI)
{  
    F_Fluid[CELLI] = vector(F_Fluid_Collect[0][myid][CELLI], F_Fluid_Collect[1][myid][CELLI], 0.0);
}
MPI_Barrier(MPI_COMM_WORLD);
//Info << "CELLS[131054] F_Fluid: " << "F_Fluid.x = " << F_Fluid[131054].x() << "  " << "F_Fluid.y = " << F_Fluid[131054].y() << endl;
//Info << "CELLS[131055] F_Fluid: " << "F_Fluid.x = " << F_Fluid[131055].x() << "  " << "F_Fluid.y = " << F_Fluid[131055].y() << endl;
//Info << "CELLS[131056] F_Fluid: " << "F_Fluid.x = " << F_Fluid[131056].x() << "  " << "F_Fluid.y = " << F_Fluid[131056].y() << endl;
