#define JIT_ELEMENT_SHARED_LIB
#define ASSEMBLE_HESSIAN_VIA_SYMMETRY
#include "jitbridge.h"

static JITFuncSpec_Table_FiniteElement_t * my_func_table;
#include "jitbridge_hang.h"


static void ResidualAndJacobian0(const JITElementInfo_t * eleminfo, const JITShapeInfo_t * shapeinfo,double * residuals, double *jacobian, double *mass_matrix,unsigned flag)
{
  int local_eqn, local_unknown;
  unsigned nummaster,nummaster2;
  double hang_weight,hang_weight2;
  const double * t=shapeinfo->t;
  const double * dt=shapeinfo->dt;
  const unsigned this_nodalind_p = 1;
  const unsigned this_nodalind_h = 0;
  //START: Precalculate time derivatives of the necessary data
  PYOOMPH_AQUIRE_ARRAY(double, this_d1t0BDF2_h, eleminfo->nnode_C2)
  for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
  {
    this_d1t0BDF2_h[l_shape]=0.0;
    for (unsigned tindex=0;tindex<shapeinfo->timestepper_ntstorage;tindex++)
    {
      this_d1t0BDF2_h[l_shape] += shapeinfo->timestepper_weights_dt_BDF2_degr[tindex]*eleminfo->nodal_data[l_shape][this_nodalind_h][tindex];
    }
  }
  //END: Precalculate time derivatives of the necessary data

  //START: Spatial integration loop
  for(unsigned ipt=0;ipt<shapeinfo->n_int_pt;ipt++)
  {
    my_func_table->fill_shape_buffer_for_point(ipt, &(my_func_table->shapes_required_ResJac[0]), flag);
    const double dx = shapeinfo->int_pt_weight;
    //START: Interpolate all required fields
    double this_intrp_d0t0_d1x1_p=0.0;
    double this_intrp_d0t0_d0x_p=0.0;
    double this_intrp_d0t0_d1x0_p=0.0;
    double this_intrp_d0t0_d1x1_h=0.0;
    double this_intrp_d0t0_d0x_h=0.0;
    double this_intrp_d0t0_d1x0_h=0.0;
    double this_intrp_d1t0BDF2_d0x_h=0.0;
    for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
    {
      this_intrp_d0t0_d1x1_p+= eleminfo->nodal_data[l_shape][this_nodalind_p][0] * shapeinfo->dx_shape_C2[l_shape][1];
      this_intrp_d0t0_d0x_p+= eleminfo->nodal_data[l_shape][this_nodalind_p][0] * shapeinfo->shape_C2[l_shape];
      this_intrp_d0t0_d1x0_p+= eleminfo->nodal_data[l_shape][this_nodalind_p][0] * shapeinfo->dx_shape_C2[l_shape][0];
      this_intrp_d0t0_d1x1_h+= eleminfo->nodal_data[l_shape][this_nodalind_h][0] * shapeinfo->dx_shape_C2[l_shape][1];
      this_intrp_d0t0_d0x_h+= eleminfo->nodal_data[l_shape][this_nodalind_h][0] * shapeinfo->shape_C2[l_shape];
      this_intrp_d0t0_d1x0_h+= eleminfo->nodal_data[l_shape][this_nodalind_h][0] * shapeinfo->dx_shape_C2[l_shape][0];
      this_intrp_d1t0BDF2_d0x_h+= this_d1t0BDF2_h[l_shape] * shapeinfo->shape_C2[l_shape];
    }
    //END: Interpolate all required fields


    // SUBEXPRESSIONS

 //Subexpressions // TODO: Check whether it is constant to take it out of the loop
    //Derivatives of subexpressions
    if (flag)
    {
    }
    //START: Contribution of the spaces
    double _res_contrib,_J_contrib;
    {
      double const * testfunction = shapeinfo->shape_C2;
      DX_SHAPE_FUNCTION_DECL(dx_testfunction) = shapeinfo->dx_shape_C2;
      DX_SHAPE_FUNCTION_DECL(dX_testfunction) = shapeinfo->dX_shape_C2;
      for (unsigned int l_test=0;l_test<eleminfo->nnode_C2;l_test++)
      {
        BEGIN_RESIDUAL_CONTINUOUS_SPACE(eleminfo->nodal_local_eqn[l_test][this_nodalind_h],dx*( 3.3333333333333331e-01*dx_testfunction[l_test][1]*this_intrp_d0t0_d1x1_p*pow(this_intrp_d0t0_d0x_h,3.0)+3.3333333333333331e-01*dx_testfunction[l_test][0]*this_intrp_d0t0_d1x0_p*pow(this_intrp_d0t0_d0x_h,3.0)+testfunction[l_test]*this_intrp_d1t0BDF2_d0x_h), shapeinfo->hanginfo_C2,this_nodalind_h,l_test)
          ADD_TO_RESIDUAL_CONTINUOUS_SPACE()
          BEGIN_JACOBIAN()
            for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
            {
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_p], ( 3.3333333333333331e-01*dx_testfunction[l_test][1]*shapeinfo->dx_shape_C2[l_shape][1]*pow(this_intrp_d0t0_d0x_h,3.0)+3.3333333333333331e-01*dx_testfunction[l_test][0]*shapeinfo->dx_shape_C2[l_shape][0]*pow(this_intrp_d0t0_d0x_h,3.0))*dx,shapeinfo->hanginfo_C2,this_nodalind_p,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_h], dx*( dx_testfunction[l_test][1]*this_intrp_d0t0_d1x1_p*pow(this_intrp_d0t0_d0x_h,2.0)*shapeinfo->shape_C2[l_shape]+dx_testfunction[l_test][0]*this_intrp_d0t0_d1x0_p*pow(this_intrp_d0t0_d0x_h,2.0)*shapeinfo->shape_C2[l_shape]+testfunction[l_test]*shapeinfo->timestepper_weights_dt_BDF2_degr[0]*shapeinfo->shape_C2[l_shape]),shapeinfo->hanginfo_C2,this_nodalind_h,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
                ADD_TO_MASS_MATRIX_HANG_HANG(testfunction[l_test]*dx*shapeinfo->shape_C2[l_shape])
              END_JACOBIAN_HANG()
            }
          END_JACOBIAN()
        END_RESIDUAL_CONTINUOUS_SPACE()
        BEGIN_RESIDUAL_CONTINUOUS_SPACE(eleminfo->nodal_local_eqn[l_test][this_nodalind_p],( testfunction[l_test]*( this_intrp_d0t0_d0x_p+-1.4457428321908239e-05*( pow(this_intrp_d0t0_d0x_h,3.0)-4.2187499999999993e-07)/pow(this_intrp_d0t0_d0x_h,6.0))-dx_testfunction[l_test][1]*this_intrp_d0t0_d1x1_h-dx_testfunction[l_test][0]*this_intrp_d0t0_d1x0_h)*dx, shapeinfo->hanginfo_C2,this_nodalind_p,l_test)
          ADD_TO_RESIDUAL_CONTINUOUS_SPACE()
          BEGIN_JACOBIAN()
            for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
            {
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_p], testfunction[l_test]*dx*shapeinfo->shape_C2[l_shape],shapeinfo->hanginfo_C2,this_nodalind_p,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_h], ( testfunction[l_test]*( 8.6744569931449428e-05*( pow(this_intrp_d0t0_d0x_h,3.0)-4.2187499999999993e-07)/pow(this_intrp_d0t0_d0x_h,7.0)*shapeinfo->shape_C2[l_shape]+-4.3372284965724714e-05*1.0/pow(this_intrp_d0t0_d0x_h,4.0)*shapeinfo->shape_C2[l_shape])-dx_testfunction[l_test][1]*shapeinfo->dx_shape_C2[l_shape][1]-dx_testfunction[l_test][0]*shapeinfo->dx_shape_C2[l_shape][0])*dx,shapeinfo->hanginfo_C2,this_nodalind_h,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
            }
          END_JACOBIAN()
        END_RESIDUAL_CONTINUOUS_SPACE()
      }
    }
    //END: Contribution of the spaces
  }
  //END: Spatial integration loop

}


// INITIAL CONDITION 
static double ElementalInitialConditions0(const JITElementInfo_t * eleminfo, int field_index,double *_x, double *_xlagr,double *_normal,double t,int flag,double default_val)
{
  if (field_index==0) // IC of field h
  {
    if (!flag) return 1.2500000000000000e+00*fmax(fmax( -4.0000000000000002e-01*(_x[1]*_x[1])+-4.0000000000000002e-01*(_x[0]*_x[0])-_x[0]-2.2500000000000001e-01,  -4.0000000000000002e-01*(_x[1]*_x[1])+-4.0000000000000002e-01*(_x[0]*_x[0])+_x[0]-2.2500000000000001e-01), 6.0000000000000001e-03); 
    if (flag==1) return 0.0; 
    if (flag==2) return 0.0; 
  }
  return default_val;
}

static double ElementalDirichletConditions(const JITElementInfo_t * eleminfo, int field_index,double *_x, double *_xlagr,double *_normal,double t,double default_val)
{
  const unsigned flag=0;
  return default_val;
}

// Used for Z2 error estimators
static double GeometricJacobian(const JITElementInfo_t * eleminfo, const double * _x)
{
  return 1.0000000000000000e+00;
}
// Used for elemsize_Eulerian etc
static double JacobianForElementSize(const JITElementInfo_t * eleminfo, const double * _x)
{
  return 1.0000000000000000e+00;
}


static void GetZ2Fluxes(const JITElementInfo_t * eleminfo, const JITShapeInfo_t * shapeinfo, double * Z2Flux)
{

  const unsigned this_nodalind_h = 0;
  //START: Precalculate time derivatives of the necessary data
  //END: Precalculate time derivatives of the necessary data

    //START: Interpolate all required fields
    double this_intrp_d0t0_d1x1_h=0.0;
    double this_intrp_d0t0_d1x0_h=0.0;
    for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
    {
      this_intrp_d0t0_d1x1_h+= eleminfo->nodal_data[l_shape][this_nodalind_h][0] * shapeinfo->dx_shape_C2[l_shape][1];
      this_intrp_d0t0_d1x0_h+= eleminfo->nodal_data[l_shape][this_nodalind_h][0] * shapeinfo->dx_shape_C2[l_shape][0];
    }
    //END: Interpolate all required fields

  Z2Flux[0] = this_intrp_d0t0_d1x0_h;
  Z2Flux[1] = this_intrp_d0t0_d1x1_h;
}



static void clean_up(JITFuncSpec_Table_FiniteElement_t *functable)
{
#ifndef NULL
#define PYOOMPH_NULL (void *)0
#else
#define PYOOMPH_NULL NULL
#endif
 pyoomph_tested_free(functable->fieldnames_Pos[ 0]); functable->fieldnames_Pos[0]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->fieldnames_Pos[ 1]); functable->fieldnames_Pos[1]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->fieldnames_Pos[ 2]); functable->fieldnames_Pos[2]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->fieldnames_Pos[ 3]); functable->fieldnames_Pos[3]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->fieldnames_Pos); functable->fieldnames_Pos=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->fieldnames_C2[ 0]); functable->fieldnames_C2[0]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->fieldnames_C2[ 1]); functable->fieldnames_C2[1]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->fieldnames_C2); functable->fieldnames_C2=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->ResidualAndJacobian_NoHang); functable->ResidualAndJacobian_NoHang=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->ResidualAndJacobian); functable->ResidualAndJacobian=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->ResidualAndJacobianSteady); functable->ResidualAndJacobianSteady=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->shapes_required_ResJac); functable->shapes_required_ResJac=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->shapes_required_Hessian); functable->shapes_required_Hessian=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->HessianVectorProduct); functable->HessianVectorProduct=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->missing_residual_assembly); functable->missing_residual_assembly=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->res_jac_names[0]); functable->res_jac_names[0]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->res_jac_names); functable->res_jac_names=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->temporal_error_scales); functable->temporal_error_scales=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->discontinuous_refinement_exponents); functable->discontinuous_refinement_exponents=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->IC_names[0]); functable->IC_names[0]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->IC_names); functable->IC_names=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->InitialConditionFunc); functable->InitialConditionFunc=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->Dirichlet_names[0]); functable->Dirichlet_names[0]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->Dirichlet_names[1]); functable->Dirichlet_names[1]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->Dirichlet_names[2]); functable->Dirichlet_names[2]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->Dirichlet_names[3]); functable->Dirichlet_names[3]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->Dirichlet_names[4]); functable->Dirichlet_names[4]=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->Dirichlet_names); functable->Dirichlet_names=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->Dirichlet_set); functable->Dirichlet_set=PYOOMPH_NULL; 
 pyoomph_tested_free(functable->domain_name); functable->domain_name=PYOOMPH_NULL; 
}

JIT_API void JIT_ELEMENT_init(JITFuncSpec_Table_FiniteElement_t *functable)
{
 functable->check_compiler_size(sizeof(char),1, "char");
 functable->check_compiler_size(sizeof(unsigned short),2, "unsigned short");
 functable->check_compiler_size(sizeof(unsigned int),4, "unsigned int");
 functable->check_compiler_size(sizeof(unsigned long int),8, "unsigned long int");
 functable->check_compiler_size(sizeof(unsigned long long int),8, "unsigned long long int");
 functable->check_compiler_size(sizeof(float),4, "float");
 functable->check_compiler_size(sizeof(double),8, "double");
 functable->check_compiler_size(sizeof(size_t),8, "size_t");
 functable->check_compiler_size(sizeof(struct JITElementInfo),96, "struct JITElementInfo");
 functable->check_compiler_size(sizeof(struct JITHangInfoEntry),16, "struct JITHangInfoEntry");
 functable->check_compiler_size(sizeof(struct JITHangInfo),16, "struct JITHangInfo");
 functable->check_compiler_size(sizeof(struct JITShapeInfo),568, "struct JITShapeInfo");
 functable->check_compiler_size(sizeof(struct JITFuncSpec_RequiredShapes_FiniteElement),48, "struct JITFuncSpec_RequiredShapes_FiniteElement");
 functable->check_compiler_size(sizeof(struct JITFuncSpec_Callback_Entry),32, "struct JITFuncSpec_Callback_Entry");
 functable->check_compiler_size(sizeof(struct JITFuncSpec_MultiRet_Entry),24, "struct JITFuncSpec_MultiRet_Entry");
 functable->check_compiler_size(sizeof(struct JITFuncSpec_Table_FiniteElement),1104, "struct JITFuncSpec_Table_FiniteElement");
 functable->nodal_dim=2;
 functable->lagr_dim=2;
 functable->fd_jacobian=false; 
 functable->fd_position_jacobian=false; 
 functable->with_adaptivity=true; 
 functable->debug_jacobian_epsilon = 0;
 functable->stop_on_jacobian_difference = false;
 functable->numfields_Pos=4;
 functable->fieldnames_Pos=(char **)malloc(sizeof(char*)*functable->numfields_Pos);
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,0, "coordinate_x" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,1, "coordinate_y" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,2, "lagrangian_x" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,3, "lagrangian_y" );
 functable->numfields_C2=2;
 functable->numfields_C2_bulk=2;
 functable->numfields_C2_basebulk=2;
 functable->numfields_C2_new=2;
 functable->nodal_offset_C2_basebulk =0;
 functable->buffer_offset_C2_basebulk =0;
 functable->fieldnames_C2=(char **)malloc(sizeof(char*)*functable->numfields_C2);
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C2,0, "h" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C2,1, "p" );
 functable->dominant_space=strdup("C2");
 functable->hangindex_Pos=-1; //Position always hangs on the max space
 functable->hangindex_C2TB=-1;
 functable->hangindex_C2=-1;
 functable->hangindex_C1TB=functable->numfields_C2TB_basebulk+functable->numfields_C2_basebulk;
 functable->hangindex_C1=functable->numfields_C2TB_basebulk+functable->numfields_C2_basebulk;
 functable->max_dt_order=1;
 functable->moving_nodes=false;
 functable->num_res_jacs=1;
 functable->ResidualAndJacobian_NoHang=(JITFuncSpec_ResidualAndJacobian_FiniteElement *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_ResidualAndJacobian_FiniteElement));
 functable->ResidualAndJacobian=(JITFuncSpec_ResidualAndJacobian_FiniteElement *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_ResidualAndJacobian_FiniteElement));
 functable->ResidualAndJacobianSteady=(JITFuncSpec_ResidualAndJacobian_FiniteElement *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_ResidualAndJacobian_FiniteElement));
 functable->shapes_required_ResJac=(JITFuncSpec_RequiredShapes_FiniteElement_t *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_RequiredShapes_FiniteElement_t));
 functable->shapes_required_Hessian=(JITFuncSpec_RequiredShapes_FiniteElement_t *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_RequiredShapes_FiniteElement_t));
 functable->HessianVectorProduct=(JITFuncSpec_HessianVectorProduct_FiniteElement *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_HessianVectorProduct_FiniteElement));
 functable->res_jac_names=(char**)calloc(functable->num_res_jacs,sizeof(char*));
 functable->missing_residual_assembly=(bool*)calloc(functable->num_res_jacs,sizeof(bool));
 SET_INTERNAL_FIELD_NAME(functable->res_jac_names,0, "" );
 functable->ResidualAndJacobian_NoHang[0]=&ResidualAndJacobian0;
 functable->ResidualAndJacobian[0]=&ResidualAndJacobian0;
 functable->ResidualAndJacobianSteady[0]=&ResidualAndJacobian0;
  functable->shapes_required_ResJac[0].dx_psi_C2 = true;
  functable->shapes_required_ResJac[0].psi_C2 = true;
 functable->missing_residual_assembly[0] = false;

 functable->num_Z2_flux_terms = 2;
 functable->GetZ2Fluxes=&GetZ2Fluxes;
  functable->shapes_required_Z2Fluxes.dx_psi_C2 = true;
 functable->temporal_error_scales=calloc(8,sizeof(double)); 
 functable->discontinuous_refinement_exponents=calloc(8,sizeof(double));
 functable->num_ICs=1;
 functable->IC_names=(char**)calloc(functable->num_ICs,sizeof(char*));
 functable->InitialConditionFunc=(JITFuncSpec_InitialCondition_FiniteElement*)calloc(functable->num_ICs,sizeof(JITFuncSpec_InitialCondition_FiniteElement));
 SET_INTERNAL_FIELD_NAME(functable->IC_names,0, "" );
 functable->InitialConditionFunc[0]=&ElementalInitialConditions0;
 functable->DirichletConditionFunc=&ElementalDirichletConditions;
 functable->Dirichlet_set_size=5;
 functable->Dirichlet_set=(bool *)calloc(functable->Dirichlet_set_size,sizeof(bool)); 
 functable->Dirichlet_names=(char**)calloc(functable->Dirichlet_set_size,sizeof(char*));
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,0, "" ); // nodal_coords index is 2
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,1, "coordinate_y" ); // nodal_coords index is 1
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,2, "coordinate_x" ); // nodal_coords index is 0
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,3, "h" ); // nodal_data index is 0
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,4, "p" ); // nodal_data index is 1
 functable->integration_order=0;
 functable->GeometricJacobian=&GeometricJacobian;
 functable->JacobianForElementSize=&JacobianForElementSize;
 SET_INTERNAL_NAME(functable->domain_name,"domain");
 functable->clean_up=&clean_up;
 my_func_table=functable;
}


