import argparse
import gzip
import itertools
import sys


def combine_barcodes_product(file_paths, output_path, sep=""):
    """读取多个 bc.txt 文件，进行笛卡尔积全组合并输出至 gzip 压缩文件。"""
    sequences_list = []

    # 1. 读取所有文件的序列列表，自动去除空白行与换行符
    for path in file_paths:
        try:
            with open(path, "r", encoding="utf-8") as f:
                seqs = [line.strip() for line in f if line.strip()]
                if not seqs:
                    print(f"警告: 文件 {path} 为空！", file=sys.stderr)
                sequences_list.append(seqs)
        except Exception as e:
            print(f"读取文件 {path} 失败: {e}", file=sys.stderr)
            sys.exit(1)

    # 2. 生成笛卡尔积全组合迭代器
    combo_generator = itertools.product(*sequences_list)

    # 确保输出文件后缀补全 .gz
    if not output_path.endswith(".gz"):
        output_path += ".gz"

    # 3. 实时压缩写入 gzip 文件
    count = 0
    try:
        with gzip.open(output_path, "wt", encoding="utf-8") as out_f:
            for combo in combo_generator:
                out_f.write(sep.join(combo) + "\n")
                count += 1
        print(
            f"成功处理 {len(file_paths)} 个文件，共生成 {count} 条笛卡尔积组合序列，已保存至: {output_path}"
        )
    except Exception as e:
        print(f"写入 gzip 文件失败: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="将 2 个或多个 bc.txt 文件中的序列进行笛卡尔积全组合，并输出为 .gz 压缩格式"
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="输入文件路径，例如: bc1.txt bc2.txt [bc3.txt ...]",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="必须指定输出 gzip 文件路径（如: combined_bc.txt.gz）",
    )
    parser.add_argument(
        "-s",
        "--sep",
        default="",
        help="序列间的分隔符（可选，默认为空；如设为 '-' 则输出 ATCG-GCTA）",
    )

    args = parser.parse_args()

    if len(args.inputs) < 2:
        print("错误：请传入至少 2 个 bc.txt 文件！", file=sys.stderr)
        sys.exit(1)

    combine_barcodes_product(
        args.inputs, output_path=args.output, sep=args.sep
    )