#' 智能安装并加载R包（自动检测+安装+加载）
#' @param pkgs 包名（字符向量，GitHub包用"user/repo"格式）
#' @param dependencies 是否安装依赖项（默认TRUE）
#' @param quiet 静默模式（默认FALSE显示过程信息）
#' @return 返回包含安装和加载状态的数据框
Packages_auto <- function(pkgs, 
																dependencies = TRUE, 
																quiet = FALSE) {
	
	# 初始化结果数据框
	result <- data.frame(
		Package = pkgs,
		Installed = FALSE,
		Loaded = FALSE,
		Source = NA_character_,
		stringsAsFactors = FALSE
	)
	
	# 检查必要工具包
	if (!requireNamespace("BiocManager", quietly = TRUE)) {
		if (!quiet) message("安装BiocManager...")
		install.packages("BiocManager", quiet = TRUE)
	}
	
	if (!requireNamespace("remotes", quietly = TRUE)) {
		if (!quiet) message("安装remotes...")
		install.packages("remotes", quiet = TRUE)
	}
	
	# 主处理循环
	for (i in seq_along(pkgs)) {
		pkg <- pkgs[i]
		original_pkg <- pkg  # 保存原始名称（对GitHub包很重要）
		
		# 处理GitHub包名（转换为纯包名）
		if (grepl("/", pkg)) {
			pkg <- strsplit(pkg, "/")[[1]][2]
		}
		
		# 检查是否已安装
		if (requireNamespace(pkg, quietly = TRUE)) {
			result[i, c("Installed", "Loaded", "Source")] <- c(TRUE, TRUE, "Pre-installed")
			if (!quiet) message("✔ ", pkg, " 已安装并加载")
			next
		}
		
		# 尝试从CRAN安装
		if (!quiet) message("尝试从CRAN安装: ", pkg)
		cran_success <- tryCatch({
			install.packages(pkg, dependencies = dependencies, quiet = TRUE)
			requireNamespace(pkg, quietly = TRUE)
		}, error = function(e) FALSE)
		
		if (cran_success) {
			result[i, c("Installed", "Loaded", "Source")] <- c(TRUE, TRUE, "CRAN")
			if (!quiet) message("✔ ", pkg, " 从CRAN安装成功")
			next
		}
		
		# 尝试从Bioconductor安装
		if (!quiet) message("尝试从Bioconductor安装: ", pkg)
		bioc_success <- tryCatch({
			BiocManager::install(pkg, dependencies = dependencies, ask = FALSE, quiet = TRUE)
			requireNamespace(pkg, quietly = TRUE)
		}, error = function(e) FALSE)
		
		if (bioc_success) {
			result[i, c("Installed", "Loaded", "Source")] <- c(TRUE, TRUE, "Bioconductor")
			if (!quiet) message("✔ ", pkg, " 从Bioconductor安装成功")
			next
		}
		
		# 尝试从GitHub安装（如果是user/repo格式）
		if (grepl("/", original_pkg)) {
			if (!quiet) message("尝试从GitHub安装: ", original_pkg)
			gh_success <- tryCatch({
				remotes::install_github(original_pkg, dependencies = dependencies, quiet = TRUE)
				requireNamespace(pkg, quietly = TRUE)
			}, error = function(e) FALSE)
			
			if (gh_success) {
				result[i, c("Installed", "Loaded", "Source")] <- c(TRUE, TRUE, "GitHub")
				if (!quiet) message("✔ ", pkg, " 从GitHub安装成功")
				next
			}
		}
		
		# 全部尝试失败
		if (!quiet) message("× ", pkg, " 安装失败")
	}
	
	# 打印汇总信息
	if (!quiet) {
		message("\n===== 最终状态 =====")
		print(result)
		
		failed <- result$Package[!result$Loaded]
		if (length(failed) > 0) {
			message("\n警告：以下包未成功加载：", paste(failed, collapse = ", "))
		}
	}
	
	invisible(result)
}


# pkgs <- c("dplyr", "Seurat", "Signac", "timoast/signac")
# status <- smart_load_packages(pkgs)