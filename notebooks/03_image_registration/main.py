import os
from valis import registration
from PIL import Image, ImageOps

sample = 'X6'
mode = ['register', 'transfer', 'all'][1]

# params for mode = 'register'
param = 2000
negate_dapi = False

# params for mode = 'transfer'
prot_channel = 'Channel02'

# 0. Set dirs
dapi_dir = f"/valis/images/{sample}/dapi/"
prot_dir = f"/valis/images/{sample}/prot/"
ref_slide = dapi_dir + f"{sample}_dapi.tiff"

out_dir = f"/valis/out/{sample}/"
reg_dir = "/valis/registered/"

# Regsiter
if mode == 'register' or mode == 'all':
    # 1. Make negative phenocycker DAPI
    pheno_dapi_path = prot_dir + f"Phenocycler_Channel01_{sample}_Level2.tif"
    pheno_dapi_path_reg = dapi_dir + f"Phenocycler_Channel01_{sample}_Level2_reg.tif"

    if negate_dapi:
        with Image.open(pheno_dapi_path) as img:
            negative_img = ImageOps.invert(img)
            negative_img.save(pheno_dapi_path_reg, format="TIFF")
    else:
        # os.system(f'cp {pheno_dapi_path} {pheno_dapi_path_reg}')
        with Image.open(pheno_dapi_path) as img:
            negative_img = img
            negative_img.save(pheno_dapi_path_reg, format="TIFF")


    print(f"Negative image saved as {pheno_dapi_path_reg}")

    # 2. Register the DAPIs
    registrar = registration.Valis(
        dapi_dir,
        out_dir,
        reference_img_f = ref_slide,
        max_processed_image_dim_px = param)

    rigid_registrar, non_rigid_registrar, error_df = registrar.register()

    registrar.warp_and_save_slides(reg_dir, pyramid=False)

# Transfer
if mode == 'transfer' or mode == 'all':
    # 3. Transform all prot chanells
    registrar_file = out_dir + "dapi/data/dapi_registrar.pickle"
    print(registrar_file)

    new_image_dir = prot_dir + f"Phenocycler_{prot_channel}_{sample}_Level2.tif"
    reg_slide_dir = reg_dir + f"Phenocycler_{prot_channel}_{sample}_Level2.tiff"

    registrar = registration.load_registrar(registrar_file)
    print(registrar.name_dict)
    slide_obj = registrar.get_slide(f'Phenocycler_Channel01_{sample}_Level2_reg')

    slide_obj.warp_and_save_slide(
        src_f=new_image_dir,
        dst_f=reg_slide_dir,
        pyramid=False,
        )

# End
print(f'process done for {sample} chennel {prot_channel}')
registration.kill_jvm()
