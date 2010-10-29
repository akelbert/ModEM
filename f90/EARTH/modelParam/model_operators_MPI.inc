
#ifdef MPI
  !Create a vector suitable for MPI messages from modelParam_t

  subroutine create_model_param_place_holder(sigma_mpi)

      type (modelParam_t), intent(in)   	        :: sigma_mpi
      integer  dimension1

       dimension1 = count(sigma_mpi%c%exists)

         if(associated(model_para_vec)) then
             deallocate(model_para_vec)
         end if
       allocate(model_para_vec(dimension1))

  end subroutine create_model_param_place_holder


   subroutine get_model_para_values(sigma_mpi)

      type (modelParam_t), intent(in)   	        :: sigma_mpi

    call getParamValues_modelParam(sigma_mpi, model_para_vec)

   end subroutine get_model_para_values


   subroutine set_model_para_values(sigma_mpi)

      type (modelParam_t), intent(inout)   	        :: sigma_mpi

    call fillParamValues_modelParam(sigma_mpi, model_para_vec)

   end subroutine set_model_para_values

#endif
