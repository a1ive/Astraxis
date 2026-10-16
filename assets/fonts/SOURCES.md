# 字体来源

| 文件 | 来源 | 许可 |
|---|---|---|
| `Jost-300-Light.ttf`、`Jost-400-Book.ttf` | Jost* 3.5（indestructible type*，Owen Earl），发行包 <https://github.com/indestructible-type/Jost/releases/download/3.5/Jost.zip> 里的 `fonts/ttf/` 下的同名文件，未修改。SHA-256：Light `88027c9e1ff51b10373b950c962d2430b38523a76496a6d8cdc65ccd99408d6e`，Book `60de951651870fd2dbbd099a96f9321f183171cd7d73047905787f37d2ec2a13` | SIL Open Font License 1.1，Copyright 2019 indestructible type*；许可证全文见同目录 `OFL.txt`（发行包里的 `LICENSE.md`） |
| `DidactGothic-Regular.ttf` | Didact Gothic（Daniel Johnson），Google Fonts 仓库 <https://github.com/google/fonts/raw/4be3f642b64e428908c337fd95a1ef79c3a6c595/ofl/didactgothic/DidactGothic-Regular.ttf>，未修改。SHA-256 `b98f9e091b6337a86f4e0f69c94c31905e1e287782853f472e176b4e28d3f81f` | SIL Open Font License 1.1，Copyright 2009–2011 Daniel Johnson；许可证全文见同目录 `DidactGothic-OFL.txt`（仓库同目录的 `OFL.txt`） |

用途（`src/app/captions.cpp`）：切换场景时的标题卡和事件展签的标题用 Light，展签正文用 Book。Jost 没有希腊字母（只有 Ω、π 等几个符号），整个希腊字母区（U+0370–03FF）改用 Didact Gothic，合并进两种 Jost，按 x 高度缩放到与 Jost 一致。界面和标签仍用 Dear ImGui 自带的字体。
