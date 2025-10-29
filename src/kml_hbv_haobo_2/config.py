from src.config.database import GENOTYPE_REFERENCE, DATABASE


def get_thread_dict(max_threads: int) -> dict:
    """获线程数字典, 最大/高/低 线程数"""
    return {
        'high': max(1, max_threads // 2),
        'low': max(1, max_threads // 8),
        'max': max_threads
    }


def get_custom_params(
        max_threads: int,
        fb_para_num: int = 0,
) -> dict:
    """
    获取自定义参数, 包括 freebayes 线程数参数
    :param max_threads:     最大线程数
    :param fb_para_num:  freebayes 线程数参数
    :return:                自定义参数字典
    """
    if fb_para_num == 0:
        fb_para_num = max_threads
    return {
        'freebayes_threads': max(1, max_threads // fb_para_num)
    }


def database_update_ref(genotype: str) -> dict:
    """获取数据库基因型参考文件路径"""
    database = DATABASE.copy()
    database['ref'] = GENOTYPE_REFERENCE[genotype]
    return database
