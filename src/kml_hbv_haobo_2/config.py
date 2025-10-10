def get_thread_dict(max_threads: int) -> dict:
    """获线程数字典, 最大/高/低 线程数"""
    return {
        'high': max(1, max_threads // 2),
        'low': max(1, max_threads // 8),
        'max': max_threads
    }


def get_custom_params(
        max_threads: int,
        freebayes_para_num: int = 0,
) -> dict:
    """
    获取自定义参数, 包括 freebayes 线程数参数
    :param max_threads:     最大线程数
    :param freebayes_para_num:  freebayes 线程数参数
    :return:                自定义参数字典
    """
    if freebayes_para_num == 0:
        freebayes_para_num = max_threads
    return {
        'freebayes_threads': max(1, max_threads // freebayes_para_num)
    }
