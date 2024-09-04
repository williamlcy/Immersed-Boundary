std::vector<scalar> P_Interpolation;

std::vector<std::vector<double>> pLocal_Collect;
pLocal_Collect.resize(mpi_num);
for(int k=0;k<mpi_num;k++)
{
        pLocal_Collect[k].resize(cellNum[k]);
}
{

std::vector<double> pLocal(local_cellnum);
for(int j=0;j<local_cellnum;j++)
{
	pLocal[j]=p[j];
}
MPI_Barrier(MPI_COMM_WORLD);
{

	if(myid!=0)
	{
		MPI_Send(pLocal.data(), local_cellnum, MPI_DOUBLE, 0, 555, MPI_COMM_WORLD);
	}
	else
	{

		pLocal_Collect[0]=pLocal;
		for(int k=1;k<mpi_num;k++)
		{
		 MPI_Recv(pLocal_Collect[k].data(), cellNum[k], MPI_DOUBLE, k, 555, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(int kk=0;kk<mpi_num;kk++)
	{
		if(myid==0)
		{
			for(int k=1;k<mpi_num;k++)
			{
			   MPI_Send(pLocal_Collect[kk].data(), cellNum[kk], MPI_DOUBLE, k, 666+kk, MPI_COMM_WORLD);
			}
		}
		else
		{
			MPI_Recv(pLocal_Collect[kk].data(), cellNum[kk], MPI_DOUBLE, 0, 666+kk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	

	MPI_Barrier(MPI_COMM_WORLD);
}
}
MPI_Barrier(MPI_COMM_WORLD);

for(label I=0; I<npt; I++)
{
    scalar P_Interpolation_Temp = 0.0; 
    label Nsupport = meshSupport[I].size();
    for(label J=0; J<Nsupport; J++)
    {
        label N = meshSupport[I][J];
        label Np = meshSupportProcessor[I][J];
        P_Interpolation_Temp += pLocal_Collect[Np][N]*Dirac[I][J]*h*h;
    }
    P_Interpolation.push_back(P_Interpolation_Temp);
}

vector pressureForce(0.0,0.0,0.0);
for(label I=0; I<npt; I++)
{
    pressureForce += P_Interpolation[I]*Lmarks_normal[I]*Lmarks_Ds;
}

Info << "Pressure Force by IBM: " << "F_x = " << pressureForce.x() << "  " << "F_y = " << pressureForce.y() << endl;
if(myid==0)
{
getchar();
}
MPI_Barrier(MPI_COMM_WORLD);
