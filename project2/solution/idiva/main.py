import argparse

from idiva.utils.clf_related import create_df, get_clf


def main(args: any):
    """
    Reads a given vcf file, transforms it into a dataframe and uses a pretrained classifier to classify each variant in
    the vcf with 0 for control and 1 for "sick". It then adds the label to the input vcf and returns it (or saves it).
    """
    df = create_df(args.vcf_file_path)
    clf = get_clf(args)
    input_data = df.to_numpy()
    prediction = clf.predict(input_data)

    # todo: somehow include the predictions into the input vcf file
    return prediction


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf_file_path', type=str, help='vcf-file that contains the variants to be classified')
    parser.add_argument('--which_clf', type=str, default='',
                        help="indicates which classifier is used to make the predictions")

    flags = parser.parse_args()
    main(flags)
