for (label level=0;level<maxRefinementLevel+1;level++)
{
    // update mesh
    Info << "Updating the mesh: initial refinment" << endl;
    mesh.update();
        
    if (mesh.changing())
    {                    
        MRF.update();

        if (correctPhi)
        {
            // Calculate absolute flux
            // from the mapped surface velocity
            phi = mesh.Sf() & Uf();

            #include "correctPhi.H"

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi, U);
        }

        if (checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }
        
        lambda *= 0; 
        HFDIBDEM.recreateBodies(lambda,refineF);
    }    
}
