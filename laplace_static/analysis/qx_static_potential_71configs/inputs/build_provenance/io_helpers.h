static int read_env_int(
    const char *name,
    const int default_value,
    const int min_value,
    const int max_value
) {
    const char *text = getenv(name);

    if (text == NULL || *text == '\0') {
        return default_value;
    }

    errno = 0;
    char *end = NULL;
    const long value = strtol(text, &end, 10);

    if (errno == ERANGE ||
        end == text ||
        *end != '\0' ||
        value < min_value ||
        value > max_value) {
        fprintf(
            stderr,
            "ERROR: %s must be an integer in [%d,%d], got '%s'\n",
            name,
            min_value,
            max_value,
            text
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    return (int)value;
}

static const char *env_or_default(
    const char *name,
    const char *default_value
) {
    const char *value = getenv(name);

    if (value == NULL || *value == '\0') {
        return default_value;
    }

    return value;
}

static void read_evec_time(
    const char *fname,
    const int nc,
    const int time,
    const char *nbase,
    qcd_geometry *geo,
    qcd_spinorComponent3d *v
) {
    DistEigvecsHdf5Reader reader;
    memset(&reader, 0, sizeof(reader));

    DistEigvecsHdf5_ComplexDouble *disteigvecs = NULL;

    int ret;
    char *tmp_run_name = NULL;
    int tmp_nc = -1;
    int tmp_time = -1;
    unsigned int tmp_nv = 0;

    int spatial_lat_size[3] = {
        geo->L[1],
        geo->L[2],
        geo->L[3]
    };

    int spatial_proclat_size[3] = {
        geo->L[1] / geo->lL[1],
        geo->L[2] / geo->lL[2],
        geo->L[3] / geo->lL[3]
    };

    int spatial_proclat_pos[3] = {
        geo->Pos[1],
        geo->Pos[2],
        geo->Pos[3]
    };

    DistEigvecsHdf5_Metadata metadata;
    DistEigvecsHdf5_Programdata programdata;

    ret = DistEigvecsHdf5Reader_Init(&reader);

    if (ret != 0) {
        fprintf(
            stderr,
            "ERROR: DistEigvecsHdf5Reader_Init failed\n"
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    ret = DistEigvecsHdf5Reader_OpenFile(
        &reader,
        geo->comm3d,
        fname,
        &tmp_run_name,
        &tmp_nc,
        spatial_lat_size,
        spatial_proclat_size,
        &tmp_time,
        &tmp_nv,
        &metadata,
        &programdata
    );

    if (ret != 0) {
        fprintf(
            stderr,
            "ERROR: OpenFile failed for %s\n",
            fname
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (metadata.bc.type!=0 || !metadata.laplacesmear.type ||
        strcmp(metadata.laplacesmear.type,"APE") || metadata.laplacesmear.itercnt!=20 ||
        !isfinite(metadata.laplacesmear.alpha) || fabs(metadata.laplacesmear.alpha-0.5)>1e-14)
        MPI_Abort(MPI_COMM_WORLD,1);
    for(int k=0;k<3;k++)if(!isfinite(metadata.bc.theta[k]) || fabs(metadata.bc.theta[k])>1e-14)MPI_Abort(MPI_COMM_WORLD,1);

    if (tmp_nc != nc ||
        tmp_time != time ||
        tmp_nv != NV ||
        tmp_run_name == NULL ||
        strcmp(tmp_run_name, nbase) != 0) {
        fprintf(
            stderr,
            "ERROR: eigenvector metadata mismatch for %s: "
            "cfg=%d time=%d Nv=%u run_name=%s\n",
            fname,
            tmp_nc,
            tmp_time,
            tmp_nv,
            tmp_run_name ? tmp_run_name : "(null)"
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    disteigvecs = malloc(
        (size_t)geo->lV3 *
        3U *
        NV *
        sizeof(*disteigvecs)
    );

    if (disteigvecs == NULL) {
        fprintf(
            stderr,
            "ERROR: allocation of distributed eigenvectors failed\n"
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    ret = DistEigvecsHdf5Reader_ReadDistEigvecsRange(
        &reader,
        spatial_proclat_pos,
        0,
        NV,
        disteigvecs
    );

    if (ret != 0) {
        fprintf(
            stderr,
            "ERROR: ReadDistEigvecsRange failed for %s\n",
            fname
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const size_t mode_bytes =
        (size_t)geo->lV3 *
        3U *
        sizeof(*disteigvecs);

    for (int j = 0; j < NV; j++) {
        memcpy(
            &(v[j].D[0][0].re),
            &(disteigvecs[(size_t)j * geo->lV3 * 3U]),
            mode_bytes
        );
    }

    free(disteigvecs);

    ret = DistEigvecsHdf5Reader_CloseFile(&reader);

    if (ret != 0) {
        fprintf(
            stderr,
            "ERROR: CloseFile failed for %s\n",
            fname
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    DistEigvecsHdf5Reader_Free(&reader);
    free(tmp_run_name);free(metadata.laplacesmear.type);
    free(programdata.name);free(programdata.infile);
}

static void temporal_product(const qcd_gaugeField *u, int x,int y,int z,
                             int first,int end,qcd_complex_16 p[3][3])
{
    mat3_zero(p);
    for(int a=0;a<3;++a) p[a][a].re=1;
    for(int t=first;t<end;++t) {
        long v=qcd_LEXIC(t,x,y,z,u->geo->lL);
        mat3_mul(p,p,u->D[v][0]);
    }
}
