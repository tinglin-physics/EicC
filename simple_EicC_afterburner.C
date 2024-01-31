// 
// Get from https://github.com/eic/afterburner and rewrite by Ting Lin using simple root macro on 01/31/2024.
// 

struct AfterburnerCorrection {
    TVector3 boost;
    TRotation rotation;
    TLorentzVector vertex;    
};

TLorentzVector boost_vector(const TLorentzVector &part_vect, const TVector3 &result_boost) {

    long double deltaX = result_boost.X();
    long double deltaY = result_boost.Y();
    long double deltaZ = result_boost.Z();
    double deltalength2 = result_boost.Mag2();
    long double deltalength = sqrt(deltalength2);
    long double gamma = 1.0/sqrt(1.0-deltalength2);

    long double tempX = part_vect.X();
    long double tempY = part_vect.Y();
    long double tempZ = part_vect.Z();
    long double tempE = part_vect.E();
    long double nr = (deltaX*tempX+deltaY*tempY+deltaZ*tempZ)/deltalength;
    long double gfac = (gamma-1)*nr/deltalength-tempE*gamma;
    tempX+=(deltaX*gfac);
    tempY+=(deltaY*gfac);
    tempZ+=(deltaZ*gfac);
    tempE = gamma*(tempE-deltalength*nr);

    return TLorentzVector(tempX, tempY, tempZ, tempE);
}  

static void RotY(double theta, double xin, double yin, double zin, double *xout, double *yout, double *zout){

    *xout = xin*cos(theta) + zin*sin(theta);
    *yout = yin;
    *zout = zin*cos(theta) - xin*sin(theta);
}

double get_collision_width(const double widthA, const double widthB){
    return widthA * widthB / sqrt(widthA * widthA + widthB * widthB);
}

double vertex_smear(const double position, const double width, const int funtype){

    assert(width >= 0);

    double res = position;

    if (width == 0)
        return res;

    if (funtype == 1) {
        return (position - width) + 2 * gRandom->Uniform(1)*width; 
    }

    if (funtype == 2) {
        return position + gRandom->Gaus(0,width);
    }

    return res;
}

TLorentzVector move_vertex(const TLorentzVector &init_vtx) {
   
    double vertex_smear_width_x = 0;
    double vertex_smear_width_y = 0;
    double vertex_smear_width_z = 0;
    double vertex_smear_width_t = 0;

    double vertex_shift_x = 0;
    double vertex_shift_y = 0;
    double vertex_shift_z = 0;
    double vertex_shift_t = 0;

    double x = init_vtx.X() + vertex_smear(vertex_shift_x, vertex_smear_width_x, 2);
    double y = init_vtx.Y() + vertex_smear(vertex_shift_y, vertex_smear_width_y, 2);
    double z = init_vtx.Z() + vertex_smear(vertex_shift_z, vertex_smear_width_z, 2);
    double t = init_vtx.T() + vertex_smear(vertex_shift_t, vertex_smear_width_t, 2);

    return TLorentzVector(x, y, z, t);
}   

TLorentzVector generate_vertx_with_bunch_interaction(const double hadron_z, const double lepton_z, const double crossing_angle, const double sigma_hor, const double sigma_ver){

    double c_c = cos(crossing_angle/2.0);
    double s_c = sin(crossing_angle/2.0);
    double t_c = tan(crossing_angle/2.0);

    double t_int = (lepton_z - hadron_z)/(2.0*c_c);
    double z_int = (lepton_z + hadron_z)/2.0;
    double z_bunch_int = c_c*t_int;
    double x_int = z_bunch_int*t_c;

    double y_int = 0.;

    if(abs(sigma_hor) > 1e-9)
    {
        x_int += gRandom->Gaus(0,sigma_hor);
    }

    if(abs(sigma_ver) > 1e-9)
    {
        y_int += gRandom->Gaus(0,sigma_ver);
    }

    // We now have the x-y-z position of the collision in the accelerator frame, but we want it in the detector frame. Rotate by 0.5*theta_c to get to accelerator frame

    double tmpVtxX, tmpVtxY, tmpVtxZ;
    tmpVtxX = tmpVtxY = tmpVtxZ = 0.;

    RotY(crossing_angle/2.0,x_int,y_int,z_int,&tmpVtxX,&tmpVtxY,&tmpVtxZ);

    TVector3 collision_center(tmpVtxX, tmpVtxY, tmpVtxZ);

    TLorentzVector result_vertex(collision_center, t_int);

    return result_vertex;
}

TVector3 smear_beam_divergence(const TVector3 &beam_dir, const double beam_divergence_hor, const double beam_divergence_ver, const double vtx_z, const double crab_hor, const double crab_ver){

    // y direction in accelerator
    static const TVector3 accelerator_plane(0, 1, 0);

    // Horizontal
    double horizontal_angle = vtx_z * crab_hor;
    horizontal_angle = horizontal_angle + gRandom->Gaus(0,beam_divergence_hor);

    TRotation x_smear_in_accelerator_plane;
    x_smear_in_accelerator_plane.Rotate(horizontal_angle,accelerator_plane);

    double vertical_angle =  -vtx_z * crab_ver;
    vertical_angle = vertical_angle + gRandom->Gaus(0,beam_divergence_ver);

    TVector3 out_accelerator_plane = accelerator_plane.Cross(beam_dir);

    TRotation y_smear_out_accelerator_plane;
    y_smear_out_accelerator_plane.Rotate(vertical_angle,out_accelerator_plane);

    TVector3 real_dir = y_smear_out_accelerator_plane * x_smear_in_accelerator_plane * beam_dir;

    return real_dir;
}

AfterburnerCorrection construct_afterburner(){

//Beam parameters
    static constexpr double nm  = 1.e-6;
    static constexpr double cm  = 10.;    

    // Crossing angle
    float crossing_angle_hor = 50e-3;          // Crossing angle in horizontal plane [rad]
    float crossing_angle_ver = 0;         // Crossing angle in vertical plane [rad]

    float hadron_beam_beta_crab_hor = 2330000.0;
    float lepton_beam_beta_crab_hor =  800000.0;

    // Beam divergence
    float hadron_beam_divergence_hor = 1.4e-3;
    float hadron_beam_divergence_ver = 2.0e-3;
    float lepton_beam_divergence_hor = 0.71e-3;
    float lepton_beam_divergence_ver = 0.61e-3;

    // Beam beta star [mm]
    float hadron_beam_beta_star_hor = 50;
    float hadron_beam_beta_star_ver = 12;
    float lepton_beam_beta_star_hor = 100;
    float lepton_beam_beta_star_ver = 40;

    // RMS emittance
    float hadron_beam_rms_emittance_hor = 100 * nm;
    float hadron_beam_rms_emittance_ver = 50 * nm;
    float lepton_beam_rms_emittance_hor = 50 * nm;
    float lepton_beam_rms_emittance_ver = 15 * nm;

    // RMS bunch length
    float hadron_beam_rms_bunch_length = 8 * cm;
    float lepton_beam_rms_bunch_length = 0.75 * cm;

//--------------------------------------------------//

// bunch interaction simulation    
    // now handle the collision vertex first, in the head-on collision frame
    // this is used as input to the Crab angle correction

    // boost-rotation from beam angles
    TLorentzVector init_vtx(0,0,0,0);

    gRandom->SetSeed(0);

    double hadron_z;
    double lepton_z;

    TLorentzVector result_vertex;

    bool use_beam_bunch_sim = 1;

    if(use_beam_bunch_sim){
	    hadron_z = gRandom->Gaus(0,hadron_beam_rms_bunch_length);
	    lepton_z = gRandom->Gaus(0,lepton_beam_rms_bunch_length);

	    double crossing_angle = crossing_angle_hor;

	    double had_bunch_rms_hor = sqrt(hadron_beam_beta_star_hor * hadron_beam_rms_emittance_hor);
	    double had_bunch_rms_ver = sqrt(hadron_beam_beta_star_ver * hadron_beam_rms_emittance_ver);
	    double lep_bunch_rms_hor = sqrt(lepton_beam_beta_star_hor * lepton_beam_rms_emittance_hor);
	    double lep_bunch_rms_ver = sqrt(lepton_beam_beta_star_ver * lepton_beam_rms_emittance_ver);
	    double sigma_hor = get_collision_width(had_bunch_rms_hor, lep_bunch_rms_hor);
	    double sigma_ver = get_collision_width(had_bunch_rms_ver, lep_bunch_rms_ver);

	    result_vertex = generate_vertx_with_bunch_interaction(hadron_z, lepton_z, crossing_angle, sigma_hor, sigma_ver);
    }else{
        result_vertex = move_vertex(init_vtx);
	hadron_z = lepton_z = result_vertex.Z();
    }

    TVector3 z_axis(0, 0, 1);
    TVector3 y_axis(0, 1, 0);

    TVector3 ideal_lepton_dir(0, 0, -1);
    TVector3 ideal_hadron_dir(0, 0, 1);
    ideal_hadron_dir.RotateY(crossing_angle_hor);
    ideal_hadron_dir.RotateX(crossing_angle_ver);

    // Calculate angular deflection
    double hadron_crab_hor = crossing_angle_hor / 2.0 / sqrt(hadron_beam_beta_crab_hor * hadron_beam_beta_star_hor);
    double hadron_crab_ver = 0;

    double lepton_crab_hor = crossing_angle_hor / 2.0 / sqrt(lepton_beam_beta_crab_hor * lepton_beam_beta_star_hor);
    double lepton_crab_ver = 0;

    //Smear hadron beam divergence
    TVector3 real_hadron_dir = smear_beam_divergence(ideal_hadron_dir, hadron_beam_divergence_hor, hadron_beam_divergence_ver, hadron_z, hadron_crab_hor, hadron_crab_ver);
    TVector3 real_lepton_dir = smear_beam_divergence(ideal_lepton_dir, lepton_beam_divergence_hor, lepton_beam_divergence_ver, lepton_z, lepton_crab_hor, lepton_crab_ver);

    TVector3 boost_axis = real_hadron_dir + real_lepton_dir;

    TVector3 result_boost = 0.5 * boost_axis;

    TVector3 beamDiffAxis = (real_hadron_dir - real_lepton_dir);

    beamDiffAxis = (1.0/beamDiffAxis.Mag())*beamDiffAxis;

    double cos_rotation_angle_to_z = beamDiffAxis.Dot(z_axis);

    TVector3 rotation_axis_1 = (real_hadron_dir - real_lepton_dir).Cross(z_axis);

    const double rotation_angle_to_z_1 = acos(cos_rotation_angle_to_z);

    TRotation result_rotation;
    result_rotation.Rotate(rotation_angle_to_z_1,rotation_axis_1);

    TVector3 beamCenterDiffAxis = (ideal_hadron_dir - ideal_lepton_dir);
    beamCenterDiffAxis = (1.0/beamCenterDiffAxis.Mag())*beamCenterDiffAxis;

    double cos_rotation_center_angle_to_z = beamCenterDiffAxis.Dot(z_axis);
    
    TVector3 rotation_axis_2 = beamCenterDiffAxis.Cross(z_axis);

    const double rotation_angle_to_z_2 = -acos(cos_rotation_center_angle_to_z);

    TRotation temp_rotation;
    temp_rotation.Rotate(rotation_angle_to_z_2,rotation_axis_2);

    TVector3 init_3vertex(result_vertex.X(), result_vertex.Y(), result_vertex.Z());

    TVector3 final_3vertex = init_3vertex;

    result_vertex.SetXYZT(final_3vertex.X(), final_3vertex.Y(), final_3vertex.Z(), result_vertex.T());

    AfterburnerCorrection result;

    result.boost = result_boost;
    result.rotation = result_rotation;
    result.vertex = result_vertex;

    return result;
}

int simple_EicC_afterburner(){

    AfterburnerCorrection ab_correction;
    ab_correction = construct_afterburner();
    Double_t angle;
    TVector3 axis;
    ab_correction.rotation.AngleAxis(angle,axis);

    TLorentzVector part_vect_ini;
    part_vect_ini.SetPxPyPzE(-8.974310634048452e-02,-1.787917727861366e+00,-2.441600654467308e+00,3.027559684589505e+00);

    cout<<"Before beam smearing: px: "<<part_vect_ini.Px()<<"; py: "<<part_vect_ini.Py()<<"; pz: "<<part_vect_ini.Pz()<<"; energy: "<<part_vect_ini.E()<<endl;

    TLorentzVector part_vect;
    part_vect = boost_vector(part_vect_ini,ab_correction.boost);
    part_vect.Rotate(angle, axis);

    cout<<"After  beam smearing: px: "<<part_vect.Px()<<"; py: "<<part_vect.Py()<<"; pz: "<<part_vect.Pz()<<"; energy: "<<part_vect.E()<<endl;

    return 1;
}	
