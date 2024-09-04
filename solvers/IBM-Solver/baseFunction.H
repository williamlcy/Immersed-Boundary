//Roma A M, Peskin C S, Berger M J.
// An adaptive version of the immersed boundary method[J].
// Journal of computational physics, 1999, 153(2): 509-534.
scalar Delta(scalar R)
{
    scalar COF;

    if(0.5 <= R && R <= 1.5)
    {
        COF = 1.0/6*(5.0 - 3.0*R - Foam::sqrt(-3.0*(1.0-R)*(1.0-R)+1.0));
    }
    else if(R <= 0.5)
    {
        COF = 1.0/3*(1 + Foam::sqrt(-3.0*R*R+1.0));
    }
    else
    {
        COF = 0;
    }
    return COF;
}