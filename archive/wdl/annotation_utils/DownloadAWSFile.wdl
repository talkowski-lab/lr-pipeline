version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow DownloadAWSFile {
    meta {
        description: [
            "This utility downloads a single file from S3 and copies it to GCS, mirroring the S3 path structure relative to a configurable base prefix."
        ]
    }

    parameter_meta {
        aws_path: "S3 URI of the file to download."
        gcs_folder: "GCS destination folder."
        base_path: "S3 base prefix to strip when constructing the destination GCS path."
        gcs_path: "GCS URI of the transferred file."
    }

    input {
        String aws_path
        String gcs_folder
        String base_path

        String utils_docker
        RuntimeAttr? runtime_attr_override
    }

    # Strip the base S3 prefix and append to GCS folder to mirror the AWS path structure
    String relative_path = sub(aws_path, "^" + base_path, "")
    String output_gcs_path = sub(gcs_folder, "/+$", "") + "/" + relative_path

    call Helpers.TransferAWSToGCS {
        input:
            aws_path = aws_path,
            output_gcs_path = output_gcs_path,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        String gcs_path = output_gcs_path
    }
}
